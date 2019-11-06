import argparse
import numpy as np
import cooler

import logging
log = logging.getLogger(__name__)
from multiprocessing import Process, Queue
import time

from krbalancing import *

from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicmatrix.lib import MatrixFileHandler


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description='Correct each matrix of the given mcool matrix with KR correction.')

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='Matrix to reduce in h5 format.',
                                metavar='matrix.h5',
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the resulting matrix. '
                                'The output is also a .h5 file. But don\'t add '
                                'the suffix.',
                                required=True)
    # parserRequired.a
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


def compute_merge(pMatrixName, pMatrixList, pQueue):

    out_queue_list = []
    for matrix in pMatrixList:
        hic = hm.hiCMatrix(pMatrixName + '::' + matrix)

        kr = kr_balancing(hic.matrix.shape[0], hic.matrix.shape[1],
                          hic.matrix.count_nonzero(), hic.matrix.indptr.astype(np.int64, copy=False),
                          hic.matrix.indices.astype(np.int64, copy=False), hic.matrix.data.astype(np.float64, copy=False))
        kr.computeKR()
        correction_factors = kr.get_normalisation_vector(False).todense()
        hic.setCorrectionFactors(correction_factors)
        out_queue_list.append(hic)

    pQueue.put(out_queue_list)
    return


def main(args=None):

    args = parse_arguments().parse_args(args)

    threads = args.threads
    merged_matrices = [None] * threads
    matrices_list = cooler.fileops.list_coolers(args.matrix)
    if len(matrices_list) > threads:
        threads = len(matrices_list)
    all_data_collected = False
    thread_done = [False] * threads
    log.debug('matrix read, starting processing')
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
        process[i] = Process(target=compute_merge, kwargs=dict(
            pMatrixName=args.matrix,
            pMatrixList=matrices_name_list,
            pQueue=queue[i]
        )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                merged_matrices[i] = queue[i].get()

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

    merged_matrices = [item for sublist in merged_matrices for item in sublist]
    log.debug('len(merged_matrices) {}'.format(len(merged_matrices)))
    log.debug('len(matrices_list) {}'.format(len(matrices_list)))

    for i, hic_matrix in enumerate(merged_matrices):
        append = False
        if i > 0:
            append = True
        matrixFileHandlerOutput = MatrixFileHandler(
            pFileType='cool', pAppend=append, pFileWasH5=False)

        matrixFileHandlerOutput.set_matrix_variables(hic_matrix.matrix, hic_matrix.cut_intervals, hic_matrix.nan_bins,
                                                     hic_matrix.correction_factors, hic_matrix.distance_counts)
        matrixFileHandlerOutput.save(args.outFileName + '::' + matrices_list[i], pSymmetric=True, pApplyCorrection=False)
