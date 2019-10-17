import argparse
import os
import gzip
import shutil
from multiprocessing import Process, Queue
import time

import logging
log = logging.getLogger(__name__)

import cooler
from sparse_neighbors_search import MinHash
from sparse_neighbors_search import MinHashSpectralClustering
from sparse_neighbors_search import MinHashClustering

from sklearn.cluster import KMeans

from hicmatrix import HiCMatrix as hm
from schicexplorer.utilities import opener
from hicmatrix.lib import MatrixFileHandler

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import skfuzzy
from scipy.spatial.distance import cdist
from sklearn.cluster import KMeans
from multiprocessing import Process, Queue
import scipy.sparse

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
                           choices=['spectral', 'kmeans', 'fuzzy-cmeans', 'mst'],
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
        _matrix = hic_ma.matrix

        if neighborhood_matrix is None:
            neighborhood_matrix = csr_matrix((pXDimension, _matrix.shape[0] * _matrix.shape[1]), dtype=np.float)

        instances, features = _matrix.nonzero()

        instances *= _matrix.shape[1]
        instances += features
        features = None
        neighborhood_matrix[pIndex+i, instances] = _matrix.data
    
    pQueue.put(neighborhood_matrix)

def fcm(data, n_clusters=1, n_init=30, m=2, max_iter=300, tol=1e-16):
    # Copied from https://codereview.stackexchange.com/questions/188455/fuzzy-c-means-in-python?rq=1
    min_cost = np.inf
    for iter_init in range(n_init):

        # Randomly initialize centers
        centers = data[np.random.choice(
            data.shape[0], size=n_clusters, replace=False
            ), :]

        # Compute initial distances
        # Zeros are replaced by eps to avoid division issues
        dist = np.fmax(
            cdist(centers, data, metric='sqeuclidean'),
            np.finfo(np.float64).eps
        )

        for iter1 in range(max_iter):

            # Compute memberships       
            u = (1 / dist) ** (1 / (m-1))
            um = (u / u.sum(axis=0))**m

            # Recompute centers
            prev_centers = centers
            centers = um.dot(data) / um.sum(axis=1)[:, None]

            dist = cdist(centers, data, metric='sqeuclidean')

            if np.linalg.norm(centers - prev_centers) < tol:
                break

        # Compute cost
        cost = np.sum(um * dist)
        if cost < min_cost:
            min_cost = cost
            min_centers = centers
            mem = um.argmax(axis=0)

    return min_centers, mem
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

    if args.clusterMethod == 'spectral':
        log.debug('spectral clustering')
        minHashSpectralClustering = MinHashSpectralClustering(n_clusters=args.numberOfClusters, number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
                                                            shingle_size=4, fast=args.fastModeMinHash, n_neighbors=neighborhood_matrix.shape[0])
        log.debug('spectral clustering fit predict')

        labels_clustering = minHashSpectralClustering.fit_predict(neighborhood_matrix)
        log.debug('create label matrix assoziation')
    elif args.clusterMethod == 'kmeans':
        log.debug('kmeans clustering')
        kmeans_object = KMeans(n_clusters=args.numberOfClusters, random_state=0, n_jobs=args.threads)
        minHash_object = MinHash(number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
                                                            shingle_size=4, fast=args.fastModeMinHash, n_neighbors=neighborhood_matrix.shape[0])
        # (n_clusters=args.numberOfClusters, number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
                                                            # shingle_size=4, fast=args.fastModeMinHash, n_neighbors=args.numberOfNeighbors)
        minHashClustering = MinHashClustering(minHashObject=minHash_object, clusteringObject=kmeans_object)
        log.debug('kmeans clustering fit predict')

        labels_clustering = minHashClustering.fit_predict(neighborhood_matrix)
    elif args.clusterMethod == 'fuzzy-cmeans':
        log.debug('start fuzzy cmeans')
        
        minHash_object = MinHash(number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
                                                            shingle_size=4, fast=args.fastModeMinHash, n_neighbors=neighborhood_matrix.shape[0])

        minHash_object.fit(neighborhood_matrix)
        precomputed_graph = minHash_object.kneighbors_graph(mode='distance')
        log.debug('minhash done, start fuzzy clustering')       
        # cntr, u, u0, d, jm, p, fpc = skfuzzy.cluster.cmeans(
        # precomputed_graph, args.numberOfClusters, 2, error=0.005, maxiter=300, init=None)
        # log.debug('cntr {} '.format(cntr))
        min_centers, mem = fcm(data=precomputed_graph.toarray(), n_clusters=args.numberOfClusters)
        log.debug('min_centers {}'.format(min_centers))
        log.debug('min_centers {} {}'.format(len(min_centers), len(min_centers[0])))
        log.debug('precomputed_graph {}'.format(precomputed_graph.shape))
        labels_clustering = min_centers.T

        log.debug('mem {}'.format(mem))
    elif args.clusterMethod == 'mst':
        minHash_object = MinHash(number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
                                                            shingle_size=4, fast=args.fastModeMinHash, n_neighbors=neighborhood_matrix.shape[0])

        minHash_object.fit(neighborhood_matrix)
        precomputed_graph = minHash_object.kneighbors_graph(mode='distance')

        import networkx 
        import matplotlib.pyplot as plt
        graph = networkx.convert_matrix.from_numpy_array(precomputed_graph.toarray())
        T = networkx.minimum_spanning_tree(graph)

        pos = networkx.spring_layout(T, iterations=100)
        # plt.subplot(211)
        # nx.draw(H, pos, with_labels=False, node_size=10)

        networkx.draw(T, pos, node_size=2)
        plt.savefig('mst.png', dpi=300)
        plt.close()

        networkx.draw_random(T, node_size=2)
        plt.savefig('random.png', dpi=300)
        plt.close()
        log.debug("plot II")

        networkx.draw_circular(T, node_size=2)
        plt.savefig('circular.png', dpi=300)
        plt.close()
        log.debug("plot III")

        networkx.draw_spectral(T, node_size=2)
        plt.savefig('spectral.png', dpi=300)
        plt.close()
        log.debug("plot IV")




    matrices_cluster = list(zip(matrices_list, labels_clustering))
    np.savetxt(args.outFileName, matrices_cluster, fmt="%s")
