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


def compute_correction(pMatrixName, pMatrixList, pQueue):

    out_queue_list = []

    index_datatype = np.int64
    valid_matrix_list = []
    try:
        for i, matrix in enumerate(pMatrixList):

            pixels_chromosome, shape, _ = load_matrix(pMatrixName + '::' + matrix, pChromosomes, False, None)

            if max_shape < shape[0]:
                max_shape = shape[0]

            _matrix = [None, None, None]
            if 'bin1_id' in pixels_chromosome.columns and 'bin2_id' in pixels_chromosome.columns and 'count' in pixels_chromosome.columns:
                _matrix[0] = pixels_chromosome['bin1_id'].values
                _matrix[1] = pixels_chromosome['bin2_id'].values
                _matrix[2] = pixels_chromosome['count'].values

                # matrix = csr_matrix()
                matrix = csr_matrix((data, (instances, features)), (pXDimension, max_shape * max_shape), dtype=np.float)
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

            # hic.setCorrectionFactors(correction_factors)
            # out_queue_list.append(hic)

    except Exception as exp:
        print('hff')

    # for matrix in pMatrixList:
    #     hic = hm.hiCMatrix(pMatrixName + '::' + matrix)

    #     kr = kr_balancing(hic.matrix.shape[0], hic.matrix.shape[1],
    #                       hic.matrix.count_nonzero(), hic.matrix.indptr.astype(np.int64, copy=False),
    #                       hic.matrix.indices.astype(np.int64, copy=False), hic.matrix.data.astype(np.float64, copy=False))
    #     kr.computeKR()
    #     correction_factors = kr.get_normalisation_vector(False).todense()
    #     hic.setCorrectionFactors(correction_factors)
    #     out_queue_list.append(hic)

    pQueue.put(out_queue_list)
    return


def main(args=None):

    args = parse_arguments().parse_args(args)

    threads = args.threads
    matrixFileHandler_list = [None] * threads
    matrices_list = cell_name_list(args.matrix)
    if len(matrices_list) < threads:
        threads = len(matrices_list)
    all_data_collected = False
    thread_done = [False] * threads
    length_index = [None] * threads
    length_index[0] = 0
    matricesPerThread = len(matrices_list) // threads
    queue = [None] * threads
    process = [None] * threads
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
            pQueue=queue[i]
        )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                matrixFileHandler_list[i] = queue[i].get()

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

    # merged_matrices = [item for sublist in merged_matrices for item in sublist]
    # matrixFileHandlerObjects_list = []
    # for i, hic_matrix in enumerate(merged_matrices):
    #     matrixFileHandlerOutput = MatrixFileHandler(pMatrixFile=matrices_list[i],
    #                                                 pFileType='cool', pFileWasH5=False)

    #     matrixFileHandlerOutput.set_matrix_variables(hic_matrix.matrix, hic_matrix.cut_intervals, hic_matrix.nan_bins,
    #                                                  hic_matrix.correction_factors, hic_matrix.distance_counts)
    #     # matrixFileHandlerOutput.save(args.outFileName + '::' + matrices_list[i], pSymmetric=True, pApplyCorrection=False)
    #     matrixFileHandlerObjects_list.append(matrixFileHandlerOutput)
    # matrixFileHandler = MatrixFileHandler(pFileType='scool')
    # matrixFileHandler.matrixFile.coolObjectsList = matrixFileHandlerObjects_list
    # matrixFileHandler.save(args.outFileName, pSymmetric=True, pApplyCorrection=False)

    matrix_file_handler_object_list = [item for sublist in matrixFileHandler_list for item in sublist]

    matrixFileHandler = MatrixFileHandler(pFileType='scool')
    matrixFileHandler.matrixFile.coolObjectsList = matrix_file_handler_object_list
    matrixFileHandler.save(args.outFileName, pSymmetric=True, pApplyCorrection=False)
