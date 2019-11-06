
import argparse
import os
from multiprocessing import Process, Queue
import time

import logging
log = logging.getLogger(__name__)

import cooler

from hicexplorer.hicAdjustMatrix import adjustMatrix
from schicexplorer._version import __version__
from hicmatrix.lib import MatrixFileHandler


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        add_help=False
    )

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='The matrix to adjust in the mcool format.',
                                required=True)
    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the adjusted matrix.',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')
    parserMutuallyExclusive = parser.add_mutually_exclusive_group()
    parserMutuallyExclusive.add_argument('--chromosomes', '-c',
                                         nargs='+',
                                         help='List of chromosomes to keep / remove')
    parserMutuallyExclusive.add_argument('--regions', '-r',
                                         help='BED file which stores a list of regions to keep / remove')
    parserMutuallyExclusive.add_argument('--maskBadRegions', '-mbr',
                                         help='Bad regions are identified and masked.')
    parserOpt.add_argument('--action',
                           help='Keep, remove or mask the list of specified chromosomes / regions ',
                           default='keep',
                           choices=['keep', 'remove', 'mask']
                           )
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def compute_adjust_matrix(pMatrixName, pMatricesList, pArgs, pQueue):
    hicmatrices_adjusted_objects = []
    for i, matrix in enumerate(pMatricesList):
        pArgs.matrix = pMatrixName + '::' + matrix
        hic_ma_adjusted = adjustMatrix(pArgs)
        hicmatrices_adjusted_objects.append(hic_ma_adjusted)

    pQueue.put(hicmatrices_adjusted_objects)
    return


def main(args=None):
    # args_string
    args = parse_arguments().parse_args(args)
    hicmatrix_adjusted_objects = []
    matrices_name = args.matrix
    threads = args.threads
    matrices_list = cooler.fileops.list_coolers(matrices_name)
    if threads > len(matrices_list):
        threads = len(matrices_list)

    all_data_collected = False
    thread_done = [False] * threads
    hicmatrix_adjusted_objects_threads = [None] * threads
    matricesPerThread = len(matrices_list) // threads
    queue = [None] * threads
    process = [None] * threads
    for i in range(threads):

        if i < threads - 1:
            matrices_name_list = matrices_list[i * matricesPerThread:(i + 1) * matricesPerThread]
        else:
            matrices_name_list = matrices_list[i * matricesPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=compute_adjust_matrix, kwargs=dict(
            pMatrixName=matrices_name,
            pMatricesList=matrices_name_list,
            pArgs=args,
            pQueue=queue[i]
        )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                hicmatrix_adjusted_objects_threads[i] = queue[i].get()
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

    # TODO: implement this!
    # hicmatrix_adjusted_objects = [item for sublist in hicmatrix_adjusted_objects_threads for item in sublist]

    # for i, matrixFileHandler in enumerate(matrixFileHandlerList):
    #     matrixFileHandlerOutput = MatrixFileHandler(pFileType='cool', pAppend=pAppend, pEnforceInteger=False, pFileWasH5=False, pHic2CoolVersion=None)

    #     matrixFileHandlerOutput.set_matrix_variables(_matrix, cut_intervals, nan_bins,
    #                                                  correction_factors, distance_counts)
    #     matrixFileHandler.save(args.outFileName + '::' + matrices_list[i], pSymmetric=True, pApplyCorrection=False)
