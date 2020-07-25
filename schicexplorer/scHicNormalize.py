import argparse
from multiprocessing import Process, Queue
import time
import logging
log = logging.getLogger(__name__)

import cooler
import numpy as np
import pandas as pd

from hicmatrix import HiCMatrix
from hicmatrix.lib import MatrixFileHandler

from schicexplorer._version import __version__
from schicexplorer.utilities import cell_name_list


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description=''
    )

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--matrix', '-m',
                                help='The single cell Hi-C interaction matrices to investigate for QC. Needs to be in scool format',
                                metavar='scool scHi-C matrix',
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name of the normalized scool matrix.',
                                required=True)

    parserRequired.add_argument('--threads', '-t',
                                help='Number of threads. Using the python multiprocessing module.',
                                required=False,
                                default=4,
                                type=int)
    parserRequired.add_argument('--normalize', '-n',
                                help='Normalize to a) all matrices to the lowest read count of the given matrices, b) all to a given read coverage value or c) to a multiplicative value',
                                choices=['smallest', 'read_count', 'multiplicative'],
                                default='smallest',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--setToZeroThreshold', '-z',
                           help='Values smaller as this threshold are set to 0.',
                           required=False,
                           default=1.0,
                           type=float)
    parserOpt.add_argument('--value', '-v', default=1,
                           type=float,
                           help='This value is used to either be interpreted as the desired read_count or the multiplicative value. This depends on the value for --normalize')
    parserOpt.add_argument('--maximumRegionToConsider',
                           help='To compute the normalization factor for the normalization mode \'smallest\' and \'read_count\', consider only this genomic distance around the diagonal.',
                           required=False,
                           default=30000000,
                           type=int)
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def compute_sum(pMatrixName, pMatricesList, pMaxDistance, pQueue):
    sum_list = []
    for i, matrix in enumerate(pMatricesList):

        matrixFileHandler = MatrixFileHandler(pFileType='cool', pMatrixFile=pMatrixName + '::' + matrix, pLoadMatrixOnly=True)
        _matrix, cut_intervals, nan_bins, \
            distance_counts, correction_factors = matrixFileHandler.load()
        # try:
        instances = _matrix[0]
        features = _matrix[1]

        distances = np.absolute(instances - features)
        mask = distances <= pMaxDistance

        sum_of_matrix = _matrix[2][mask].sum()
        # instances = _matrix[0]
        # features = _matrix[1]

        # distances = np.absolute(instances - features)
        # mask = distances == 0

        # # remove the double diagonal values
        # sum_of_matrix -= _matrix[2][mask].sum()
        # except:
        # sum_list.append()

        # instances = _matrix[0]
        # features = _matrix[1]

        # distances = np.absolute(instances - features)
        # mask = distances <= max_distance
        # sparsity_length = len(_matrix[2][mask])

        # sparsity.append(sparsity_length / (shape_x * max_distance))

        sum_list.append(sum_of_matrix)
        del _matrix
    pQueue.put(sum_list)


def compute_normalize(pMatrixName, pMatricesList, pNormalizeMax, pSumOfAll, pThreshold, pMultiplicative, pQueue):

    pixelList = []
    log.debug('multiplicative {}'.format(pMultiplicative))

    for i, matrix in enumerate(pMatricesList):

        matrixFileHandler = MatrixFileHandler(pFileType='cool', pMatrixFile=pMatrixName + '::' + matrix, pLoadMatrixOnly=True)
        _matrix, cut_intervals, nan_bins, \
            distance_counts, correction_factors = matrixFileHandler.load()

        data = np.array(_matrix[2]).astype(np.float32)
        instances = np.array(_matrix[0])
        features = np.array(_matrix[1])

        mask = np.isnan(data)
        data[mask] = 0
        mask = np.isinf(data)
        data[mask] = 0

        if pMultiplicative is None:
            adjust_factor = pSumOfAll[i] / pNormalizeMax
        else:
            adjust_factor = pMultiplicative

        # log.debug('pSumOfAll[i] {}'.format(pSumOfAll[i]))
        # log.debug('pNormalizeMax {}'.format(pNormalizeMax))
        # log.debug('adjust_factor {}'.format(adjust_factor))

        if pMultiplicative is None:
            data /= adjust_factor
        else:
            data *= adjust_factor

        mask = np.isnan(data)
        data[mask] = 0

        mask = np.isinf(data)
        data[mask] = 0

        mask = data < pThreshold
        data[mask] = 0

        mask = data == 0
        instances = instances[~mask]
        features = features[~mask]
        data = data[~mask]

        pixels = pd.DataFrame({'bin1_id': instances, 'bin2_id': features, 'count': data})
        # log.debug('dtyoes: {}'.format(pixels.dtypes))
        pixelList.append(pixels)

    pQueue.put(pixelList)


def main(args=None):

    args = parse_arguments().parse_args(args)

    matrices_name = args.matrix
    matrices_list = cell_name_list(matrices_name)
    # log.debug('size of matrices_list: {}'.format(len(matrices_list)))

    threads = args.threads

    # get bin size
    cooler_obj = cooler.Cooler(matrices_name + '::' + matrices_list[0])
    bin_size = cooler_obj.binsize

    sum_list_threads = [None] * args.threads
    process = [None] * args.threads
    queue = [None] * args.threads

    # all_threads_done = False
    thread_done = [False] * args.threads
    # count_output = 0
    # count_call_of_read_input = 0
    # computed_pairs = 0
    matricesPerThread = len(matrices_list) // threads

    for i in range(args.threads):
        if i < threads - 1:
            matrices_name_list = matrices_list[i * matricesPerThread:(i + 1) * matricesPerThread]
        else:
            matrices_name_list = matrices_list[i * matricesPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=compute_sum, kwargs=dict(
            pMatrixName=matrices_name,
            pMatricesList=matrices_name_list,
            pMaxDistance=args.maximumRegionToConsider // bin_size,
            pQueue=queue[i]
        )
        )
        process[i].start()

    all_data_collected = False
    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                sum_list_threads[i] = queue[i].get()
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

    sum_of_all = [item for sublist in sum_list_threads for item in sublist]
    sum_of_all = np.array(sum_of_all)
    foo = sum_of_all[sum_of_all < 100000]
    log.debug('len(foo_ {}'.format(len(foo)))
    # log.debug('size of sum_all: {}'.format(len(sum_of_all)))
    argmin = np.argmin(sum_of_all)
    if args.normalize == 'smallest':
        normalizeMax = sum_of_all[argmin]
        multiplicative = None
    elif args.normalize == 'read_count':
        normalizeMax = args.value
        multiplicative = None

    else:
        normalizeMax = None
        multiplicative = args.value
        log.debug('multiplicative')

    log.debug('sum_of_all[:10] {}'.format(sum_of_all[:10]))
    log.debug('argmin {}'.format(argmin))
    log.debug('argminSum {}'.format(normalizeMax))

    matricesPerThread = len(matrices_list) // threads

    pixelsListThreads = [None] * args.threads
    process = [None] * args.threads
    # all_data_processed = False
    queue = [None] * args.threads

    thread_done = [False] * args.threads
    # args.threads = 1
    for i in range(args.threads):
        if i < args.threads - 1:
            matrices_name_list = matrices_list[i * matricesPerThread:(i + 1) * matricesPerThread]
            sum_of_all_list = sum_of_all[i * matricesPerThread:(i + 1) * matricesPerThread]
        else:
            matrices_name_list = matrices_list[i * matricesPerThread:]
            sum_of_all_list = sum_of_all[i * matricesPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=compute_normalize, kwargs=dict(
            pMatrixName=matrices_name,
            pMatricesList=matrices_name_list,
            pNormalizeMax=normalizeMax,
            pSumOfAll=sum_of_all_list,
            pThreshold=args.setToZeroThreshold,
            pMultiplicative=multiplicative,
            pQueue=queue[i]
        )
        )
        process[i].start()
    all_data_collected = False
    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                pixelsListThreads[i] = queue[i].get()

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

    pixelsList = [item for sublist in pixelsListThreads for item in sublist]

    cooler_obj_external = cooler.Cooler(matrices_name + '::' + matrices_list[0])
    bins = cooler_obj_external.bins()[:]

    matrixFileHandler = MatrixFileHandler(pFileType='scool')
    matrixFileHandler.matrixFile.coolObjectsList = None
    matrixFileHandler.matrixFile.bins = bins
    matrixFileHandler.matrixFile.pixel_list = pixelsList
    matrixFileHandler.matrixFile.name_list = matrices_list
    matrixFileHandler.save(args.outFileName, pSymmetric=True, pApplyCorrection=False)
