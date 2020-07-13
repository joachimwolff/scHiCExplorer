import argparse
import os
import time
from multiprocessing import Process, Queue


import numpy as np
from scipy.sparse import csr_matrix, vstack, save_npz, lil_matrix, dok_matrix

import logging
log = logging.getLogger(__name__)

import cooler
from sklearn.cluster import SpectralClustering, KMeans

from sparse_neighbors_search import MinHash
from sparse_neighbors_search import MinHashClustering


from hicmatrix import HiCMatrix as hm

from schicexplorer._version import __version__
from schicexplorer.utilities import cell_name_list, create_csr_matrix_all_cells

import time

def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description='scHicClusterMinHash uses kmeans or spectral clustering to associate each cell to a cluster and therefore to its cell cycle. '
        'The clustering is applied on dimension reduced data based on an approximate kNN search with the local sensitive hashing technique MinHash. This approach reduces the number of dimensions from samples * (number of bins)^2 to samples * samples. '
        'Please consider also the other clustering and dimension reduction approaches of the scHicExplorer suite. They can give you better results, '
        'can be faster or less memory demanding.'

    )

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--matrix', '-m',
                                help='The single cell Hi-C interaction matrices to cluster. Needs to be in scool format',
                                metavar='scool scHi-C matrix',
                                required=True)
    parserRequired.add_argument('--numberOfClusters', '-c',
                                help='Number of to be computed clusters',
                                required=False,
                                default=12,
                                type=int)
    parserRequired.add_argument('--clusterMethod', '-cm',
                                help='Algorithm to cluster the Hi-C matrices',
                                choices=['spectral', 'kmeans'],
                                default='spectral')
    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the resulting clusters',
                                required=True,
                                default='clusters.txt')
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--exactModeMinHash', '-em',
                           help='This option increases the runtime significantly, from a few minutes to half an hour or longer. If set, the number of hash collisions is only used for candidate set creation and the euclidean distance is considered too.',
                           action='store_false')
    parserOpt.add_argument('--saveIntermediateRawMatrix', '-sm',
                           help='This option activates the save of the intermediate raw scHi-C matrix.',
                           required=False)
    parserOpt.add_argument('--numberOfHashFunctions', '-nh',
                           help='Number of to be used hash functions for minHash',
                           required=False,
                           default=800,
                           type=int)
    parserOpt.add_argument('--numberOfNearestNeighbors', '-k',
                           help='Number of to be used computed nearest neighbors for the knn graph.',
                           required=False,
                           default=100,
                           type=int)
    parserOpt.add_argument('--shareOfMatrixToBeTransferred', '-s',
                           help='Which share of rows shall be transferred from Python to C++ at once. Values between 0 and 1, the more are transferred at once, the larger the memory usage is. The less rows are transferred, the slower the computation is.',
                           required=False,
                           default=0.25,
                           type=float)
    parserOpt.add_argument('--chromosomes',
                           help='List of to be plotted chromosomes',
                           nargs='+')

    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module.',
                           required=False,
                           default=4,
                           type=int)
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


# def open_and_store_matrix(pMatrixName, pMatricesList, pIndex, pXDimension, pChromosomes, pQueue):
#     neighborhood_matrix = None
#     time_load = 0.0
#     time_all = 0.0
#     time_csr_create = 0.0
#     time_add = 0.0
#     features = []
#     data = []
#     features_length = []
#     max_shape = 0
#     index_datatype = np.int64
#     for i, matrix in enumerate(pMatricesList):
#         time_start_load = time.time()
#         time_start_all = time.time() 

#         if pChromosomes is not None and len(pChromosomes) == 1:
#             hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName + '::' + matrix, pChrnameList=pChromosomes, pNoIntervalTree=True, pUpperTriangleOnly=True, pLoadMatrixOnly=True, pRestoreMaskedBins=False)
#         else:
#             if not pChromosomes:
#                 hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName + '::' + matrix, pNoIntervalTree=True, pUpperTriangleOnly=True, pLoadMatrixOnly=True, pRestoreMaskedBins=False)
#             else:
#                 hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName + '::' + matrix, pNoIntervalTree=False, pUpperTriangleOnly=True, pLoadMatrixOnly=True, pRestoreMaskedBins=False)
#             if pChromosomes:
#                 hic_ma.keepOnlyTheseChr(pChromosomes)
#         _matrix = hic_ma.matrix
#         time_load += time.time() - time_start_load
#         time_csr_create_start = time.time()

#         time_csr_create += time.time() - time_csr_create_start 
#         time_add_start = time.time()
#         if max_shape < _matrix[3]:
#             max_shape = _matrix[3]
        

#         _matrix[0] = _matrix[0].astype(index_datatype)
#         _matrix[1] = _matrix[1].astype(index_datatype)

#         _matrix[0] *= np.int64(_matrix[3]) # matrix[0] are the instance ids, matrix[3] is the shape
#         _matrix[0] += _matrix[1] # matrix[3] is the shape, matrix[1] are the feature ids
#         features.extend(_matrix[0])
#         _matrix[1] = None

#         data.extend(_matrix[2])
#         features_length.append(len(_matrix[2]))
#         time_add += time.time() - time_add_start
#         del _matrix
#         time_all += time.time() - time_start_all
    

#     time_start_tocsr = time.time()
#     instances = []
#     for i, instance_length in enumerate(features_length):
#         instances.extend([pIndex + i] * instance_length)
    

#     neighborhood_matrix = csr_matrix((data, (instances, features)),(pXDimension, max_shape * max_shape), dtype=np.float)
#     log.debug('time_all {}, time_csr {}, time_add {} time_load {} time_tocsr {}'.format(time_all, time_csr_create, time_add, time_load, time.time() - time_start_tocsr))
#     pQueue.put(neighborhood_matrix)


def main(args=None):

    args = parse_arguments().parse_args(args)

    # matrices_name = args.matrix
    # threads = args.threads
    # matrices_list = cell_name_list(matrices_name)
    # neighborhood_matrix = None
    # neighborhood_matrix_threads = [None] * threads

    # all_data_collected = False
    # thread_done = [False] * threads
    # log.debug('matrix read, starting processing')
    # length_index = [None] * threads
    # length_index[0] = 0
    # matricesPerThread = len(matrices_list) // threads
    # queue = [None] * threads
    # process = [None] * threads
    # for i in range(threads):

    #     if i < threads - 1:
    #         matrices_name_list = matrices_list[i * matricesPerThread:(i + 1) * matricesPerThread]
    #         length_index[i + 1] = length_index[i] + len(matrices_name_list)
    #     else:
    #         matrices_name_list = matrices_list[i * matricesPerThread:]

    #     queue[i] = Queue()
    #     process[i] = Process(target=open_and_store_matrix, kwargs=dict(
    #         pMatrixName=matrices_name,
    #         pMatricesList=matrices_name_list,
    #         pIndex=length_index[i],
    #         pXDimension=len(matrices_list),
    #         pChromosomes=args.chromosomes,
    #         pQueue=queue[i]
    #     )
    #     )

    #     process[i].start()

    # while not all_data_collected:
    #     for i in range(threads):
    #         if queue[i] is not None and not queue[i].empty():
    #             csr_matrix_worker = queue[i].get()
    #             neighborhood_matrix_threads = csr_matrix_worker
    #             if neighborhood_matrix is None:
    #                 neighborhood_matrix = csr_matrix((len(matrices_list), neighborhood_matrix_threads.shape[1]))
    #                 neighborhood_matrix += neighborhood_matrix_threads
    #             else:
    #                 neighborhood_matrix += neighborhood_matrix_threads
    #             del neighborhood_matrix_threads
    #             queue[i] = None
    #             process[i].join()
    #             process[i].terminate()
    #             process[i] = None
    #             thread_done[i] = True
    #     all_data_collected = True
    #     for thread in thread_done:
    #         if not thread:
    #             all_data_collected = False
    #     time.sleep(1)

    # neighborhood_matrix = neighborhood_matrix_threads[0]
    # for i in range(1, len(neighborhood_matrix_threads)):
    #     neighborhood_matrix += neighborhood_matrix_threads[i]

    # neighborhood_matrix = neighborhood_matrix.tocsr()

    neighborhood_matrix, matrices_list = create_csr_matrix_all_cells(args.matrix, args.threads, args.chromosomes)
    if args.saveIntermediateRawMatrix:
        save_npz(args.saveIntermediateRawMatrix, neighborhood_matrix)
    if args.clusterMethod == 'spectral':
        log.debug('spectral clustering')
        spectral_object = SpectralClustering(n_clusters=args.numberOfClusters, affinity='nearest_neighbors', n_jobs=args.threads, random_state=0)
        log.debug('spectral clustering fit predict')
        minHash_object = MinHash(n_neighbors=args.numberOfNearestNeighbors, number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
                                 shingle_size=4, fast=args.exactModeMinHash, maxFeatures=int(max(neighborhood_matrix.getnnz(1))))
        minHashClustering = MinHashClustering(minHashObject=minHash_object, clusteringObject=spectral_object)

        labels_clustering = minHashClustering.fit_predict(neighborhood_matrix, pSaveMemory=args.shareOfMatrixToBeTransferred)
        log.debug('create label matrix assoziation')
    elif args.clusterMethod == 'kmeans':
        log.debug('kmeans clustering')
        kmeans_object = KMeans(n_clusters=args.numberOfClusters, random_state=0, n_jobs=args.threads, precompute_distances=True)
        minHash_object = MinHash(n_neighbors=args.numberOfNearestNeighbors, number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
                                 shingle_size=4, fast=args.exactModeMinHash, maxFeatures=int(max(neighborhood_matrix.getnnz(1))))
        minHashClustering = MinHashClustering(minHashObject=minHash_object, clusteringObject=kmeans_object)
        log.debug('kmeans clustering fit predict')

        labels_clustering = minHashClustering.fit_predict(neighborhood_matrix, pSaveMemory=args.shareOfMatrixToBeTransferred)

    matrices_cluster = list(zip(matrices_list, labels_clustering))
    np.savetxt(args.outFileName, matrices_cluster, fmt="%s")
