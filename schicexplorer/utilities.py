import gzip
import cooler
import os
import time
from multiprocessing import Process, Queue
import numpy as np
import pandas as pd
import logging
log = logging.getLogger(__name__)
from hicmatrix import HiCMatrix as hm
from schicexplorer._version import __version__

from scipy.sparse import csr_matrix


def opener(filename):
    """
    Determines if a file is compressed or not
    """

    f = open(filename, 'rb')

    if f.read(2) == b'\x1f\x8b':
        f.seek(0)
        return gzip.GzipFile(fileobj=f)
    else:
        f.seek(0)
        return f


def cell_name_list(pScoolUri):

    if cooler.fileops.is_scool_file(pScoolUri):
        matrices_list = cooler.fileops.list_scool_cells(pScoolUri)
        return matrices_list
    else:
        try:
            matrices_list = cooler.fileops.list_coolers(pScoolUri)

            # old and non-standard scool format stored all cells in root
            # no '/' in matrices_list and no '/cell/*'
            log.warning('Please update the scool file to the new file format standard!')
            if '/' not in matrices_list:
                return matrices_list
            # new standard scool format, all cells are stored under '/cell/'

            raise Exception('Wrong data format. Please use a scool file.')

        except Exception:
            raise Exception('Wrong data format. Please use a scool file.')
            exit(1)


def open_and_store_matrix(pMatrixName, pMatricesList, pIndex, pXDimension, pChromosomes, pIntraChromosomalContactsOnly, pChromosomeIndices, pQueue=None):
    neighborhood_matrix = None
    features = []
    data = []
    features_length = []
    max_shape = 0
    index_datatype = np.int64
    valid_matrix_list = []
    for i, matrix in enumerate(pMatricesList):

        cooler_obj = cooler.Cooler(pMatrixName + '::' + matrix)
        shape = cooler_obj.shape
        chromosome_dataframes_list = []

        if pChromosomes is None:
            pChromosomes = cooler_obj.chromnames
        for chromosome in pChromosomes:

            pixels_chromosome = cooler_obj.pixels().fetch(chromosome)
            # get pixels from one chromosome, but only intra chromosomal contacts

            # per defintion the bin1_id should belong only to the fetched chromosome, therefore only bin2_id needs to be cleaned
            if pIntraChromosomalContactsOnly:
                mask_chromosome = pixels_chromosome['bin2_id'].apply(lambda x: x in pChromosomeIndices[chromosome])
                chromosome_dataframes_list.append(pixels_chromosome[mask_chromosome])
            else:
                chromosome_dataframes_list.append(pixels_chromosome)

        pixels_chromosome = pd.concat(chromosome_dataframes_list)

        if max_shape < shape[0]:
            max_shape = shape[0]

        _matrix = [None, None, None]
        if 'bin1_id' in pixels_chromosome.columns and 'bin2_id' in pixels_chromosome.columns and 'count' in pixels_chromosome.columns:
            _matrix[0] = pixels_chromosome['bin1_id'].values
            _matrix[1] = pixels_chromosome['bin2_id'].values
            _matrix[2] = pixels_chromosome['count'].values
            if len(_matrix[2]) == 0:
                valid_matrix_list.append(False)
                continue
            valid_matrix_list.append(True)
            _matrix[0] = _matrix[0].astype(index_datatype)
            _matrix[1] = _matrix[1].astype(index_datatype)

            _matrix[0] *= np.int64(shape[0])  # matrix[0] are the instance ids, matrix[3] is the shape
            _matrix[0] += _matrix[1]  # matrix[3] is the shape, matrix[1] are the feature ids
            features.extend(_matrix[0])
            _matrix[1] = None

            data.extend(_matrix[2])
            features_length.append(len(_matrix[2]))

    instances = []
    for i, instance_length in enumerate(features_length):
        instances.extend([pIndex + i] * instance_length)

    neighborhood_matrix = csr_matrix((data, (instances, features)), (pXDimension, max_shape * max_shape), dtype=np.float)
    if pQueue is None:
        return neighborhood_matrix, valid_matrix_list
    pQueue.put([neighborhood_matrix, valid_matrix_list])
    return


def create_csr_matrix_all_cells(pMatrixName, pThreads, pChromosomes, pOutputFolder, pRawFileName, pIntraChromosomalContactsOnly=None):

    matrices_name = pMatrixName
    threads = pThreads
    matrices_list = cell_name_list(matrices_name)
    neighborhood_matrix = None
    neighborhood_matrix_threads = [None] * threads
    valid_matrix_list_threads = [None] * threads

    all_data_collected = False
    thread_done = [False] * threads
    length_index = [None] * threads
    length_index[0] = 0
    matricesPerThread = len(matrices_list) // threads
    queue = [None] * threads
    process = [None] * threads

    chromosome_indices = None
    if pIntraChromosomalContactsOnly:
        cooler_obj = cooler.Cooler(pMatrixName + '::' + matrices_list[0])
        binsDataFrame = cooler_obj.bins()[:]
        chromosome_indices = {}
        for chromosome in cooler_obj.chromnames:
            chromosome_indices[chromosome] = np.array(binsDataFrame.index[binsDataFrame['chrom'] == chromosome].tolist())

    for i in range(threads):

        if i < threads - 1:
            matrices_name_list = matrices_list[i * matricesPerThread:(i + 1) * matricesPerThread]
            length_index[i + 1] = length_index[i] + len(matrices_name_list)
        else:
            matrices_name_list = matrices_list[i * matricesPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=open_and_store_matrix, kwargs=dict(
            pMatrixName=matrices_name,
            pMatricesList=matrices_name_list,
            pIndex=length_index[i],
            pXDimension=len(matrices_list),
            pChromosomes=pChromosomes,
            pIntraChromosomalContactsOnly=pIntraChromosomalContactsOnly,
            pChromosomeIndices=chromosome_indices,
            pQueue=queue[i]
        )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                csr_matrix_worker = queue[i].get()
                neighborhood_matrix_threads, valid_matrix_list_threads[i] = csr_matrix_worker
                if neighborhood_matrix is None:
                    neighborhood_matrix = csr_matrix((len(matrices_list), neighborhood_matrix_threads.shape[1]))
                    neighborhood_matrix += neighborhood_matrix_threads
                else:
                    neighborhood_matrix += neighborhood_matrix_threads
                del neighborhood_matrix_threads
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
    return neighborhood_matrix, matrices_list



# def create_csr_matrix_share_of_cells(pMatrixName, pMatrixNameList, pChromosomes, pStartIndex, pShareSize, pLengthIndex, pIntraChromosomalContactsOnly=None):

#     # matrices_name = pMatrixName
#     # threads = pThreads
#     # # matrices_list = cell_name_list(matrices_name)
#     # neighborhood_matrix = None
#     # neighborhood_matrix_threads = [None] * threads
#     # valid_matrix_list_threads = [None] * threads

#     # all_data_collected = False
#     # thread_done = [False] * threads
#     # length_index = [None] * threads
#     # length_index[0] = 0
#     # matricesPerThread = len(matrices_list) // threads
#     # queue = [None] * threads
#     # process = [None] * threads

#     chromosome_indices = None
#     if pIntraChromosomalContactsOnly:
#         cooler_obj = cooler.Cooler(pMatrixName + '::' + pMatrixNameList[0])
#         binsDataFrame = cooler_obj.bins()[:]
#         chromosome_indices = {}
#         for chromosome in cooler_obj.chromnames:
#             chromosome_indices[chromosome] = np.array(binsDataFrame.index[binsDataFrame['chrom'] == chromosome].tolist())

#     # for i in range(threads

#     return open_and_store_matrix(
#         pMatrixName=pMatrixName,
#         pMatricesList=pMatrixNameList,
#         pIndex=pLengthIndex,
#         pXDimension=len(matrices_list),
#         pChromosomes=pChromosomes,
#         pIntraChromosomalContactsOnly=pIntraChromosomalContactsOnly,
#         pChromosomeIndices=chromosome_indices,
#         pQueue=None)
      

#     # while not all_data_collected:
#     #     for i in range(threads):
#     #         if queue[i] is not None and not queue[i].empty():
#     #             csr_matrix_worker = queue[i].get()
#     #             neighborhood_matrix_threads, valid_matrix_list_threads[i] = csr_matrix_worker
#     #             if neighborhood_matrix is None:
#     #                 neighborhood_matrix = csr_matrix((len(matrices_list), neighborhood_matrix_threads.shape[1]))
#     #                 neighborhood_matrix += neighborhood_matrix_threads
#     #             else:
#     #                 neighborhood_matrix += neighborhood_matrix_threads
#     #             del neighborhood_matrix_threads
#     #             queue[i] = None
#     #             process[i].join()
#     #             process[i].terminate()
#     #             process[i] = None
#     #             thread_done[i] = True
#     #     all_data_collected = True
#     #     for thread in thread_done:
#     #         if not thread:
#     #             all_data_collected = False
#     #     time.sleep(1)
#     return neighborhood_matrix, matrices_list