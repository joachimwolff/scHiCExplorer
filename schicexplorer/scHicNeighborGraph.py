# read all matrices
## get non-zeros and flatten it (x*length) + y
# make number of instacnes * dim**2 csr matrix

import argparse
import os
import gzip
import shutil
from multiprocessing import Process, Queue
import time
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
# warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
# 

import logging
logging.basicConfig(level=logging.DEBUG)
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('cooler').setLevel(logging.WARNING)
logging.getLogger('hicmatrix').setLevel(logging.WARNING)


log = logging.getLogger(__name__)

import cooler
from sparse_neighbors_search import MinHash
from sparse_neighbors_search import MinHashDBSCAN
from sparse_neighbors_search import MinHashSpectralClustering
from sparse_neighbors_search import MinHashClustering

from sklearn.cluster import KMeans

from hicmatrix import HiCMatrix as hm
from schicexplorer.utilities import opener
from hicmatrix.lib import MatrixFileHandler

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.decomposition import TruncatedSVD

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Process, Queue
import scipy.sparse

import networkx as nx

def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        add_help=False
    )

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--matrix', '-m',
                                help='The single cell Hi-C interaction matrices to cluster. Needs to be in mcool format',
                                metavar='mcool scHi-C matrix',
                                required=True)
    parserRequired.add_argument('--numberOfClusters', '-c',
                           help='Number of to be computed clusters',
                           required=False,
                           default=12,
                           type=int)
    parserRequired.add_argument('--numberOfNeighbors', '-n',
                           help='Number of neighbors of clustering',
                           required=False,
                           default=10,
                           type=int)
    parserRequired.add_argument('--fastModeMinHash', '-f',
                                help='If set to, only the number of hash collisions is considered for nearest neighbors search.'
                                'When not set, the number of collisions is only used for candidate set creation and the euclidean distance is considered too.',
                                action='store_true')
    parserRequired.add_argument('--numberOfHashFunctions', '-h',
                           help='Number of to be used hash functions for minHash',
                           required=False,
                           default=2000,
                           type=int)
    parserRequired.add_argument('--chromosomes',
                           help='List of to be plotted chromosomes',
                           nargs='+')
    parserRequired.add_argument('--clusterMethod', '-cm',
                           help='Algorithm to cluster the Hi-C matrices',
                           choices=['spectral', 'kmeans'],
                           default='spectral')
    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the resulting clusters',
                                required=True,
                                default='clusters.txt')
    parserRequired.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module.',
                           required=False,
                           default=4,
                           type=int)

    return parser



def open_and_store_matrix(pMatrixName, pMatricesList, pIndex, pXDimension, pChromosomes, pQueue):
    neighborhood_matrix = None
    for i, matrix in enumerate(pMatricesList):
        if pChromosomes is not None and len(pChromosomes) == 1:
            hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName+ '::' +matrix, pChrnameList=pChromosomes)
        else:
            hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName+ '::' +matrix)
            if pChromosomes:
                hic_ma.keepOnlyTheseChr(pChromosomes)
        # matrixFileHandlerInput = MatrixFileHandler(pFileType='cool', pMatrixFile=pMatrixName+ '::' +matrix)
        # _matrix, _, _, _, _ = matrixFileHandlerInput.load()
        _matrix = hic_ma.matrix

        if neighborhood_matrix is None:
            neighborhood_matrix = csr_matrix((pXDimension, _matrix.shape[0] * _matrix.shape[1]), dtype=np.float)

        instances, features = _matrix.nonzero()

        instances *= _matrix.shape[1]
        instances += features
        features = None
        neighborhood_matrix[pIndex+i, instances] = _matrix.data
        # if i % 20 == 0:
        #     log.debug('pIndex + i {} {}'.format(pIndex, i))
    
    pQueue.put(neighborhood_matrix)

def main(args=None):

    args = parse_arguments().parse_args(args)

    create_or_load_matrix = False

    # if args.createMatrix:
        
    matrices_name = args.matrix
    threads = args.threads
    matrices_list = cooler.fileops.list_coolers(matrices_name)
    neighborhood_matrix = None

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
        process[i] = Process(target=open_and_store_matrix, kwargs=dict(
                            pMatrixName = matrices_name,
                            pMatricesList= matrices_name_list, 
                            pIndex = length_index[i], 
                            pXDimension=len(matrices_list),
                            pChromosomes=args.chromosomes,
                            pQueue=queue[i]
            )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                csr_matrix_worker = queue[i].get()
                if neighborhood_matrix is None:
                    neighborhood_matrix = csr_matrix_worker
                    # log.debug('returned first csr i {}'.format(i))
                else:
                    neighborhood_matrix += csr_matrix_worker
                    # log.debug('adding csr i {}'.format(i))

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
    # scipy.sparse.save_npz(args.matrix + '_binary.npz', neighborhood_matrix)
    # else:
    #     log.debug('read npz')
    #     neighborhood_matrix = scipy.sparse.load_npz(args.matrixNpz)
    #     matrices_list = cooler.fileops.list_coolers(args.matrix)

    minhash = MinHash(number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
                        shingle_size=4, fast=args.fastModeMinHash, n_neighbors=args.numberOfNeighbors)
    
    minhash.fit(neighborhood_matrix)
    graph = minhash.kneighbors()
    log.debug('graph {}'.format(graph))
    # np.savetxt(args.outFileName, graph, fmt="%s")
    G=nx.Graph()
    for i, node in enumerate(graph[1]):
        for j, edge in enumerate(node):
            G.add_edge(i, edge, weight=graph[0][i][j])
    
    nx.draw(G, node_size=10)
    plt.savefig(args.outFileName, dpi=300)
    plt.close()
    nx.draw_random(G, node_size=10)
    plt.savefig('random_'+args.outFileName, dpi=300)
    plt.close()
    nx.draw_circular(G, node_size=10)
    plt.savefig('circular_'+args.outFileName, dpi=300)
    plt.close()
    nx.draw_spectral(G, node_size=10)
    plt.savefig('spectral_'+args.outFileName, dpi=300)
    plt.close()
    # if args.clusterMethod == 'spectral':
    #     log.debug('spectral clustering')
    #     minHashSpectralClustering = MinHashSpectralClustering(n_clusters=args.numberOfClusters, number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
    #                                                         shingle_size=4, fast=args.fastModeMinHash, n_neighbors=args.numberOfNeighbors)
    #     log.debug('spectral clustering fit predict')

    #     labels_clustering = minHashSpectralClustering.fit_predict(neighborhood_matrix)
    #     log.debug('create label matrix assoziation')
    # elif args.clusterMethod == 'kmeans':
    #     log.debug('spectral clustering')
    #     kmeans_object = KMeans(n_clusters=args.numberOfClusters, random_state=0, n_jobs=args.threads)
    #     minHash_object = MinHash(number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
    #                                                         shingle_size=4, fast=args.fastModeMinHash, n_neighbors=args.numberOfNeighbors)
    #     # (n_clusters=args.numberOfClusters, number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
    #                                                         # shingle_size=4, fast=args.fastModeMinHash, n_neighbors=args.numberOfNeighbors)
    #     minHashClustering = MinHashClustering(minHashObject=minHash_object, clusteringObject=kmeans_object)
    #     log.debug('spectral clustering fit predict')

    #     labels_clustering = minHashClustering.fit_predict(neighborhood_matrix)

    # matrices_cluster = list(zip(matrices_list, labels_clustering))
    # np.savetxt(args.outFileName, matrices_cluster, fmt="%s")
