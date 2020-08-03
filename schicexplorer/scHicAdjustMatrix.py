
import argparse
import os
from multiprocessing import Process, Queue
import time

import logging
log = logging.getLogger(__name__)

import cooler
import numpy as np
import gc
from hicmatrix import HiCMatrix as hm
import cooler
import pandas as pd
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


def compute_adjust_matrix(pMatrixName, pMatricesList, pArgs, pListIds, pInvertedMap, pInvertedLogic, pQueue):
    hicmatrices_adjusted_objects = []
    pixels_list = []
    keep_matrices = []

    for i, matrix in enumerate(pMatricesList):
        if i % 1 == 0:
            log.debug('processed {} out of {}'.format(i, len(pMatricesList)))
        try:
            cooler_obj = cooler.Cooler(pMatrixName + '::' + matrix)

            pixels = cooler_obj.pixels()[:]
            indices = pixels['bin1_id'].apply(lambda x: x in pListIds)
            if pInvertedLogic:
                indices = ~indices
            pixels = pixels[indices].reset_index(drop=True)

            indices = pixels['bin2_id'].apply(lambda x: x in pListIds)
            if pInvertedLogic:
                indices = ~indices
            pixels = pixels[indices].reset_index(drop=True)

            for key, value in pInvertedMap.items():

                pixels['bin1_id'].replace(to_replace=key, value=value, inplace=True)
                pixels['bin2_id'].replace(to_replace=key, value=value, inplace=True)

            pixels_list.append(pixels)
            keep_matrices.append(True)

        except Exception as e:
            keep_matrices.append(False)
            log.debug('exception: {}'.format(e))
            log.debug('pixels {}'.format(pixels[:5]))
            continue

    pQueue.put([pixels_list, keep_matrices])
    return


def main(args=None):
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
    if threads > len(matrices_list):
        threads = len(matrices_list)

    # load bin ids only once
    cooler_obj_external = cooler.Cooler(matrices_name + '::' + matrices_list[0])
    bins = cooler_obj_external.bins()[:]

    # apply the inverted operation if the number of values is less
    # the idea is that for
    # indices = pixels['bin1_id'].apply(lambda x: x in pListIds)
    # the search time is less if the list pListIds is shorter
    # therefore the drop must be inverted too
    apply_inverted = False
    if args.action == 'keep':
        list_ids = bins.index[bins['chrom'].apply(lambda x: x in args.chromosomes)].tolist()
        list_inverted_logic_ids = bins.index[bins['chrom'].apply(lambda x: x not in args.chromosomes)].tolist()

        bins_new = bins[bins['chrom'].apply(lambda x: x in args.chromosomes)].reset_index()

    else:
        list_ids = bins.index[bins['chrom'].apply(lambda x: x not in args.chromosomes)].tolist()
        list_inverted_logic_ids = bins.index[bins['chrom'].apply(lambda x: x in args.chromosomes)].tolist()
        bins_new = bins[bins['chrom'].apply(lambda x: x not in args.chromosomes)].reset_index()

    if len(list_inverted_logic_ids) < len(list_ids):
        apply_inverted = True
        list_ids = list_inverted_logic_ids

    dict_values = bins_new['index'].to_dict()
    inv_map = {}
    for k, v in dict_values.items():
        if k == v:
            continue
        inv_map[v] = k
    bins_new.drop(['index'], axis=1, inplace=True)

    all_data_collected = False
    thread_done = [False] * threads
    bins_thread = [None] * threads
    pixels_thread = [None] * threads
    keep_thread = [None] * threads

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
            pListIds=list_ids,
            pInvertedMap=inv_map,
            pInvertedLogic=apply_inverted,
            pQueue=queue[i]
        )
        )

        process[i].start()
    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                pixels_thread[i], keep_thread[i] = queue[i].get()

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

    pixels_list = [item for sublist in pixels_thread for item in sublist]
    keep_list = [item for sublist in keep_thread for item in sublist]

    matrices_list = np.array(matrices_list)
    mask = np.array(keep_list)
    matrices_list = matrices_list[mask]

    matrixFileHandler = MatrixFileHandler(pFileType='scool')
    matrixFileHandler.matrixFile.bins = bins_new
    matrixFileHandler.matrixFile.pixel_list = pixels_list
    matrixFileHandler.matrixFile.name_list = matrices_list

    matrixFileHandler.save(args.outFileName, pSymmetric=True, pApplyCorrection=False)
    broken_count = input_count_matrices - np.sum(np.array(keep_list))
    print('Out of {} matrices, {} were removed because they were broken.'.format(input_count_matrices, broken_count))
