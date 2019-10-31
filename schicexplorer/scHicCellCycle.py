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
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.neighbors import NearestNeighbors


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
import networkx as nx
import matplotlib.pyplot as plt

from scipy.cluster import hierarchy
from scipy.spatial import distance

from collections import defaultdict

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


def main(args=None):

    args = parse_arguments().parse_args(args)

    create_or_load_matrix = False

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
                else:
                    neighborhood_matrix += csr_matrix_worker

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

    minHash_object = MinHash(number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
                                                        shingle_size=4, fast=True, n_neighbors=args.numberOfNeighbors)

    minHash_object.fit(neighborhood_matrix)
    precomputed_graph = minHash_object.kneighbors_graph(mode='distance')

    # np.save('precomputed_graph.npy', precomputed_graph.toarray())
    graph = nx.convert_matrix.from_numpy_array(precomputed_graph.toarray())
    G = nx.minimum_spanning_tree(graph, weight='weight')
    highest_degree = sorted(list(G.degree()), key = lambda x: int(x[1]), reverse=True)[:1]

    from copy import deepcopy
    circles = []
    circles_weight = []

    # for i in range(len(precomputed_graph)):
    i = highest_degree[0][0]
    edge_list_bfs = list(nx.dfs_edges(graph, source=i, depth_limit=precomputed_graph.shape[0]))
    edge_long_list = []
    key = i
    # first_key = i
    for edge in edge_list_bfs:
        if edge[0] == key:
            edge_long_list.append([edge[0], edge[1]])
        else:
            new_append = None
            for edge_list in edge_long_list:
                if edge[0] == edge_list[-1] and edge[1] not in edge_list:
                    new_append = deepcopy(edge_list)
                    break
    #             elif edge[0] == edge_list[-1] and edge[1] in edge_list:
    #                 edge_list.append(edge[1])
    #                 circles.append(edge_list)
            if new_append is not None:
                new_append.append(edge[1])
                edge_long_list.append(new_append)

    for found_paths in edge_long_list:
        if found_paths[-1] in graph[key]:
            # found_paths.append(key)
            weight = 0
            i = 0
            while i < len(found_paths) - 1:
                weight += graph.edges[found_paths[i], found_paths[i + 1]]['weight']
                i += 1
            weight += graph.edges[found_paths[-1], key]['weight']
            weight /= (len(found_paths)+1)
            if len(found_paths) > precomputed_graph.shape[0] * 0.66:
                circles_weight.append(weight)
                circles.append(found_paths)


    # max_path = None
    # max_length = 0
    # for path in circles:
    #     if max_length < len(path):
    #         max_length = len(path)
    #         max_path = path
    max_path = circles[np.argmin(circles_weight)]
    
    matrices_list_cycle = list(np.array(matrices_list)[max_path])
    cluster_cycle = [0] * len(matrices_list_cycle)
    non_cycle = list(range(precomputed_graph.shape[0]))
    for node in max_path:
        if node in non_cycle:
            non_cycle.remove(node)

    precomputed_graph = precomputed_graph.toarray()
    nbrs = NearestNeighbors(n_neighbors=10, algorithm='ball_tree').fit(precomputed_graph)
    non_inserted = []
    for node in non_cycle:
        distances, indices = nbrs.kneighbors([precomputed_graph[node]])
        count = 0
        for neighbor in indices[0]:
            if neighbor in max_path:
    #             max_path.index(neighbor)
                max_path.insert(max_path.index(neighbor)+1, node)
                break
            else:
                count += 1
        if count == 5:
            non_inserted.append(node)

    if len(non_inserted) > 0:
        log.debug('non_inserted elements')
        matrices_list_non_cycle = list(np.array(matrices_list)[non_inserted])
        cluster_non_cycle = [1] * len(matrices_list_non_cycle)
        matrices_list_cycle.extend(matrices_list_non_cycle)
        cluster_cycle.extend(cluster_non_cycle)
        matrices_list = matrices_list_cycle
        labels_clustering = cluster_cycle
    else:
        log.debug('full list')

        matrices_list = list(np.array(matrices_list)[max_path])
        labels_clustering = [0] * len(max_path)
    
    

    matrices_cluster = list(zip(matrices_list, labels_clustering))
    np.savetxt(args.outFileName, matrices_cluster, fmt="%s")
