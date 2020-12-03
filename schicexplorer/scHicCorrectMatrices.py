import argparse
import numpy as np
import cooler

import logging
log = logging.getLogger(__name__)
logging.getLogger('hicmatrix').setLevel(logging.ERROR)

from multiprocessing import Process, Queue
import time

from krbalancing import *

from hicmatrix import HiCMatrix as hm
# from hicexplorer._version import __version__
from hicmatrix.lib import MatrixFileHandler
from schicexplorer._version import __version__
from schicexplorer.utilities import cell_name_list
from scipy.sparse import csr_matrix
from schicexplorer.utilities import load_matrix


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description='Correct each matrix of the given scool matrix with KR correction.')

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='Matrix to reduce in h5 format.',
                                metavar='matrix.h5',
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the resulting matrix, please add the scool prefix.',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module.',
                           required=False,
                           default=4,
                           type=int)
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def compute_correction(pMatrixName, pMatrixList, pCutIntervals, pQueue):

    out_queue_list = []

    print('len(pMatrixList): ' + str(len(pMatrixList)))
    try:
        for i, matrix in enumerate(pMatrixList):

            pixels, shape, _ = load_matrix(pMatrixName + '::' + matrix, None, False, None)

            # _matrix = [None, None, None]
            if 'bin1_id' in pixels.columns and 'bin2_id' in pixels.columns and 'count' in pixels.columns:
                instances = pixels['bin1_id'].values
                features = pixels['bin2_id'].values
                data = pixels['count'].values

                matrix = csr_matrix((data, (instances, features)), (shape[0], shape[1]), dtype=np.float)
            else:
                continue

            kr = kr_balancing(shape[0], shape[1],
                              matrix.count_nonzero(), matrix.indptr.astype(np.int64, copy=False),
                              matrix.indices.astype(np.int64, copy=False), matrix.data.astype(np.float64, copy=False))
            kr.computeKR()
            correction_factors = kr.get_normalisation_vector(False).todense()

            matrixFileHandlerOutput = MatrixFileHandler(pFileType='cool', pMatrixFile=matrix)

            matrixFileHandlerOutput.set_matrix_variables(matrix,
                                                         pCutIntervals,
                                                         None,
                                                         correction_factors,
                                                         None)

            out_queue_list.append(matrixFileHandlerOutput)
            print('DOne i: ' + str(i))
    except Exception as exp:
        print('Exception: ' + str(exp))
        log.debug('Exception! {}'.format(str(exp)))
        pQueue.put(str(exp))
        return

    pQueue.put(out_queue_list)
    return


def main(args=None):

    args = parse_arguments().parse_args(args)

    threads = args.threads
    matrixFileHandler_list = [None] * threads
    matrices_list = cell_name_list(args.matrix)
    if len(matrices_list) < threads:
        threads = len(matrices_list)

    matrixFileHandlerInput = MatrixFileHandler(pFileType='cool', pMatrixFile=args.matrix + "::" + matrices_list[0])

    _matrix, cut_intervals_all, nan_bins, \
        distance_counts, correction_factors = matrixFileHandlerInput.load()

    all_data_collected = False
    thread_done = [False] * threads
    length_index = [None] * threads
    length_index[0] = 0
    matricesPerThread = len(matrices_list) // threads
    queue = [None] * threads
    process = [None] * threads
    print('Threads: ' + str(threads))
    for i in range(threads):

        if i < threads - 1:
            matrices_name_list = matrices_list[i * matricesPerThread:(i + 1) * matricesPerThread]
            length_index[i + 1] = length_index[i] + len(matrices_name_list)
        else:
            matrices_name_list = matrices_list[i * matricesPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=compute_correction, kwargs=dict(
            pMatrixName=args.matrix,
            pMatrixList=matrices_name_list,
            pCutIntervals=cut_intervals_all,
            pQueue=queue[i]
        )
        )

        process[i].start()

    fail_flag = False
    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                matrixFileHandler_list[i] = queue[i].get()
                # csr_matrix_worker = queue[i].get()
                if isinstance(matrixFileHandler_list[i], str):
                    log.error('{}'.format(matrixFileHandler_list[i]))
                    fail_flag = True
                queue[i] = None
                process[i].join()
                process[i].terminate()
                process[i] = None
                thread_done[i] = True
        all_data_collected = True
        for thread in thread_done:
            if not thread:
                all_data_collected = False
        time.sleep(1)

    if fail_flag:
        exit(1)
    matrix_file_handler_object_list = [item for sublist in matrixFileHandler_list for item in sublist]

    matrixFileHandler = MatrixFileHandler(pFileType='scool')
    matrixFileHandler.matrixFile.coolObjectsList = matrix_file_handler_object_list
    matrixFileHandler.save(args.outFileName, pSymmetric=True, pApplyCorrection=False)
