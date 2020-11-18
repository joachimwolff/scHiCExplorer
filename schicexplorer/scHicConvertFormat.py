
import argparse
import os
from multiprocessing import Process, Queue
import time
import errno
import logging
log = logging.getLogger(__name__)

import cooler
import numpy as np
import gc
from hicmatrix import HiCMatrix as hm
import pandas as pd
from schicexplorer._version import __version__
from hicmatrix.lib import MatrixFileHandler

from copy import deepcopy

from schicexplorer.utilities import cell_name_list


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description='scHicConvertFormat is a tool to convert a scool matrix to other single-cell Hi-C formats. '
                    'So far only the structure of scHiCluster is supported: https://github.com/zhoujt1994/scHiCluster'
    )

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='The matrix to adjust in the scool format.',
                                required=True)
    parserRequired.add_argument('--outputFolder', '-of',
                                help='Folder name to save the files',
                                required=True,
                                default='.',
                                type=str)
    parserRequired.add_argument('--outputCellNameFile', '-oc',
                                help='File name to save the cell names and their location',
                                required=False,
                                default='cellNameFile.txt',
                                type=str)
    parserRequired.add_argument('--outputChromosomeSize', '-os',
                                help='File name to save the chromosome sizes',
                                required=False,
                                default='chromosomeSize.txt',
                                type=str)
    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--format', '-f',
                           help='The format of the output files',
                           choices=['schicluster', 'sparse-matrix-files'],
                           default='none')
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module.',
                           required=False,
                           default=4,
                           type=int)
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def convert_files(pMatrixName, pMatricesList, pBinsDataFrame, pOutputFolder, pFormat, pQueue):

    if pFormat == 'schicluster':
        return convert_to_schicluster(pMatrixName, pMatricesList, pBinsDataFrame, pOutputFolder, pQueue)
    elif pFormat == 'sparse-matrix-files':
        return convert_to_txt_csr(pMatrixName, pMatricesList, pBinsDataFrame, pOutputFolder, pQueue)
    else:
        log.error('Format not known!')
        pQueue.put([None])


def convert_to_txt_csr(pMatrixName, pMatricesList, pBinsDataFrame, pOutputFolder, pQueue):
    cell_name_array = []

    for i, matrix in enumerate(pMatricesList):

        try:
            cooler_obj = cooler.Cooler(pMatrixName + '::' + matrix)

            pixels = cooler_obj.pixels()[:]

            file_name = pOutputFolder + matrix + '.txt'
            pixels.to_csv(file_name, sep="\t", index=False, header=False)
            cell_name_array.append(pOutputFolder + matrix)

        except Exception as e:
            log.debug('exception: {}'.format(e))
            log.debug('pixels {}'.format(pixels[:5]))
            continue

    pQueue.put(cell_name_array)
    return


def convert_to_schicluster(pMatrixName, pMatricesList, pBinsDataFrame, pOutputFolder, pQueue):
    cell_name_array = []

    for i, matrix in enumerate(pMatricesList):

        try:
            cooler_obj = cooler.Cooler(pMatrixName + '::' + matrix)

            pixels = cooler_obj.pixels()[:]

            chromosome_indices = None
            for chromosome in cooler_obj.chromnames:
                chromosome_indices = np.array(pBinsDataFrame.index[pBinsDataFrame['chrom'] == chromosome].tolist())

                # get pixels from one chromosome, but only intra chromosomal contacts
                mask = pixels['bin1_id'].apply(lambda x: x in chromosome_indices) & pixels['bin2_id'].apply(lambda x: x in chromosome_indices)
                pixels_chromosome = pixels[mask].reset_index(drop=True)
                pixels_chromosome['bin1_id'] = pixels_chromosome['bin1_id'] - chromosome_indices[0]
                pixels_chromosome['bin2_id'] = pixels_chromosome['bin2_id'] - chromosome_indices[0]

                file_name = pOutputFolder + matrix + '_' + chromosome + '.txt'
                pixels_chromosome.to_csv(file_name, sep="\t", index=False, header=False)
            cell_name_array.append(pOutputFolder + matrix)

        except Exception as e:
            log.debug('exception: {}'.format(e))
            log.debug('pixels {}'.format(pixels[:5]))
            continue

    pQueue.put(cell_name_array)
    return


def main(args=None):
    args = parse_arguments().parse_args(args)
    matrices_name = args.matrix
    threads = args.threads
    matrices_list = cell_name_list(matrices_name)

    if not os.path.exists(args.outputFolder + '/cells'):
        try:
            os.makedirs(args.outputFolder + '/cells')
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    if threads > len(matrices_list):
        threads = len(matrices_list)
    # load bin ids only once
    cooler_obj = cooler.Cooler(matrices_name + '::' + matrices_list[0])
    bins = cooler_obj.bins()[:]

    all_data_collected = False
    thread_done = [False] * threads
    cell_name_array_thread = [None] * threads

    matricesPerThread = len(matrices_list) // threads
    queue = [None] * threads
    process = [None] * threads
    for i in range(threads):

        if i < threads - 1:
            matrices_name_list = matrices_list[i * matricesPerThread:(i + 1) * matricesPerThread]
        else:
            matrices_name_list = matrices_list[i * matricesPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=convert_files, kwargs=dict(
            pMatrixName=matrices_name,
            pMatricesList=matrices_name_list,
            pBinsDataFrame=bins,
            pOutputFolder=args.outputFolder,
            pFormat=args.format,
            pQueue=queue[i]
        )
        )

        process[i].start()
    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                cell_name_array_thread[i] = queue[i].get()

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

    for subset in cell_name_array_thread:
        if subset[0] is None:
            exit(1)
    cell_name_array = [item for sublist in cell_name_array_thread for item in sublist]

    # write cell names to file
    with open(args.outputCellNameFile, 'w') as file:
        for cell_name in cell_name_array:
            file.write('{}\n'.format(cell_name))
    # write chromsizes to file
    with open(args.outputChromosomeSize, 'w') as file:
        for chromosome_name, size in cooler_obj.chromsizes.items():
            file.write('{}\t{}\n'.format(chromosome_name, size))
