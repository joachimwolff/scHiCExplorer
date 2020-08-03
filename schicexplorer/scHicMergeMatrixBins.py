import argparse
import numpy as np
import cooler

import logging
log = logging.getLogger(__name__)
from multiprocessing import Process, Queue
import time

from hicmatrix import HiCMatrix as hm
from hicexplorer.hicMergeMatrixBins import running_window_merge, merge_bins
from schicexplorer._version import __version__
from hicmatrix.lib import MatrixFileHandler
from schicexplorer.utilities import cell_name_list


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description='Merges bins from a Hi-C matrix. For example, '
        'using a matrix containing 5kb bins, a matrix '
        'of 50kb bins can be derived using --numBins 10. '
    )

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='Matrix to reduce in scool format.',
                                metavar='matrix.scool',
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the resulting matrix. '
                                'The output is also a .scool file. But don\'t add '
                                'the suffix.',
                                required=True)

    parserRequired.add_argument('--numBins', '-nb',
                                help='Number of bins to merge.',
                                metavar='int',
                                type=int,
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--runningWindow',
                           help='set to merge for using a running '
                           'window of length --numBins. Must be an odd number.',
                           action='store_true')
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module.',
                           required=False,
                           default=4,
                           type=int)
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def compute_merge(pMatrixName, pMatrixList, pRunningWindow, pNumBins, pQueue):

    out_queue_list = []
    for matrix in pMatrixList:
        hic = hm.hiCMatrix(pMatrixName + '::' + matrix)

        if pRunningWindow:
            merged_matrix = running_window_merge(hic, pNumBins)
        else:
            merged_matrix = merge_bins(hic, pNumBins)

        matrixFileHandlerOutput = MatrixFileHandler(pFileType='cool', pMatrixFile=matrix, pAppend=append, pEnforceInteger=False, pFileWasH5=False)

        matrixFileHandlerOutput.set_matrix_variables(merged_matrix.matrix, merged_matrix.cut_intervals, merged_matrix.nan_bins,
                                                     merged_matrix.correction_factors, merged_matrix.distance_counts)
        out_queue_list.append(matrixFileHandlerOutput)

    pQueue.put(out_queue_list)
    return


def main(args=None):

    args = parse_arguments().parse_args(args)

    threads = args.threads
    merged_matrices = [None] * threads
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
        process[i] = Process(target=compute_merge, kwargs=dict(
            pMatrixName=args.matrix,
            pMatrixList=matrices_name_list,
            pRunningWindow=args.runningWindow,
            pNumBins=args.numBins,
            pQueue=queue[i]
        )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                log.debug('i {}'.format(i))
                log.debug('len(queue) {}'.format(len(queue)))
                log.debug('len(merged_matrices) {}'.format(len(merged_matrices)))

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

    matrixFileHandlerObjects_list = [item for sublist in merged_matrices for item in sublist]

    matrixFileHandler = MatrixFileHandler(pFileType='scool')
    matrixFileHandler.matrixFile.coolObjectsList = matrixFileHandlerObjects_list
    matrixFileHandler.save(args.outFileName, pSymmetric=True, pApplyCorrection=False)
