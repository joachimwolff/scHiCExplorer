
import argparse
import os
from multiprocessing import Process, Queue
import time

import logging
log = logging.getLogger(__name__)

import cooler
import numpy as np

from hicmatrix import HiCMatrix as hm

from schicexplorer._version import __version__
from hicmatrix.lib import MatrixFileHandler

from copy import deepcopy

from schicexplorer.utilities import cell_name_list



def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description='scHicAdjustMatrix is a tool to keep or remove a list of chromosomes of all Hi-C matrices stored in the scool file.'
    )

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='The matrix to adjust in the scool format.',
                                required=True)
    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the adjusted matrix.',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--chromosomes', '-c',
                           nargs='+',
                           help='List of chromosomes to keep / remove')
    parserOpt.add_argument('--createSubmatrix', '-cs',
                           type=int,
                           help='Keep only first n matrices and remove the rest. Good for test data creation.')
    parserOpt.add_argument('--action',
                           help='Keep, remove or mask the list of specified chromosomes / regions ',
                           default='keep',
                           choices=['keep', 'remove']
                           )
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module.',
                           required=False,
                           default=4,
                           type=int)
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def compute_adjust_matrix(pMatrixName, pMatricesList, pArgs, pQueue):
    log.debug('compute adjust matrix')
    hicmatrices_adjusted_objects = []
    keep_matrices = []
    # pArgs_local = deepcopy(pArgs)
    for i, matrix in enumerate(pMatricesList):

        hic_matrix = hm.hiCMatrix(pMatrixName + '::' + matrix)

        if pArgs.action == 'keep':
            try:
                hic_matrix.keepOnlyTheseChr(pArgs.chromosomes)
                keep_matrices.append(1)
                # hicmatrices_adjusted_objects.append(hic_matrix)
            except Exception as e:
                keep_matrices.append(0)
                log.debug('exception: {}'.format(e))
        else:
            chromosomes_list = list(hic_matrix.chrBinBoundaries)
            keep_list = []
            for chromosome in chromosomes_list:
                if chromosome in pArgs.chromosomes:
                    continue
                else:
                    keep_list.append(chromosome)
            try:
                hic_matrix.keepOnlyTheseChr(keep_list)
                keep_matrices.append(1)
            except Exception as e:
                keep_matrices.append(0)
                log.debug('exception: {}'.format(e))
        hicmatrices_adjusted_objects.append(hic_matrix)

    pQueue.put([hicmatrices_adjusted_objects, keep_matrices])
    return


def main(args=None):
    # args_string
    args = parse_arguments().parse_args(args)
    hicmatrix_adjusted_objects = []
    matrices_name = args.matrix
    threads = args.threads
    matrices_list = cell_name_list(matrices_name)
    if args.createSubmatrix is not None and args.regions is None and args.chromosomes is None:
        for matrix in matrices_list[:args.createSubmatrix]:
            cooler.fileops.cp(args.matrix + '::' + matrix, args.outFileName + '::' + matrix)
        exit(0)

    input_count_matrices = len(matrices_list)
    # log.debug('args.createSubmatrix {}, args.action {}, args.chromosomes {}'.format(args.createSubmatrix, args.action, args.chromosomes ))
    # exit()
    if threads > len(matrices_list):
        threads = len(matrices_list)

    all_data_collected = False
    thread_done = [False] * threads
    hicmatrix_adjusted_objects_threads = [None] * threads
    keep_matrices_list_threads = [None] * threads

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
    log.debug("foo")
    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                hicmatrix_adjusted_objects_threads[i], keep_matrices_list_threads[i] = queue[i].get()

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
    hicmatrix_adjusted_objects = [item for sublist in hicmatrix_adjusted_objects_threads for item in sublist]
    keep_matrices_list = [item for sublist in keep_matrices_list_threads for item in sublist]

    # matrix_names_list = []
    matrix_file_handler_object_list = []
    log.debug('length out {}'.format(len(hicmatrix_adjusted_objects)))
    for i, hic_matrix in enumerate(hicmatrix_adjusted_objects):
        if args.createSubmatrix and i > args.createSubmatrix:
            break
        append = True
        if i == 0:
            append = False

        if keep_matrices_list[i] == 0:
            continue

        matrixFileHandlerOutput = MatrixFileHandler(pFileType='cool', pAppend=append, pMatrixFile=matrices_list[i], pEnforceInteger=False, pFileWasH5=False, pHic2CoolVersion=None)

        matrixFileHandlerOutput.set_matrix_variables(hic_matrix.matrix, hic_matrix.cut_intervals, hic_matrix.nan_bins,
                                                     hic_matrix.correction_factors, hic_matrix.distance_counts)
        # matrix_names_list.append(keep_matrices_list[i])
        matrix_file_handler_object_list.append(matrixFileHandlerOutput)
        # matrixFileHandlerOutput.save(args.outFileName + '::' + matrices_list[i], pSymmetric=True, pApplyCorrection=False)

    matrixFileHandler = MatrixFileHandler(pFileType='scool')
    matrixFileHandler.matrixFile.coolObjectsList = matrix_file_handler_object_list
    matrixFileHandler.save(args.outFileName, pSymmetric=True, pApplyCorrection=False)
    broken_count = input_count_matrices - np.sum(np.array(keep_matrices_list))
    print('Out of {} matrices, {} were removed because they were broken.'.format(input_count_matrices, broken_count))
