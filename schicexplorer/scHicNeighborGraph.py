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
import networkx.algorithms.tournament as nat

def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        add_help=False
    )

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--matrix', '-m',
                                help='The single cell Hi-C interaction matrices to cluster. Needs to be in mcool format',
                                # metavar='mcool scHi-C matrix',
                                required=True)
    parserRequired.add_argument('--consensusMatrix', '-cm',
                                help='The single cell Hi-C interaction consensus matrices of the clusters. Needs to be in mcool format',
                                # metavar='mcool consensus scHi-C matrix',
                                required=True)
    # parserRequired.add_argument('--clusterFile', '-cf',
    #                             help='File containing the matrix and cluster associations.',
    #                             required=True,
    #                             default='clusters.txt')
    parserRequired.add_argument('--chromosomes',
                           help='List of to be plotted chromosomes',
                           nargs='+')
    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the resulting clusters',
                                required=True,
                                default='clusters_graph.txt')
    parserRequired.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module.',
                           required=False,
                           default=4,
                           type=int)

    return parser



def compute_consensus_sample_difference(pMatrixName, pMatricesList, pMatrixNameConsensus, pMatrixNameConsensusList, pChromosomes, pQueue):
    difference_sample_consensus_list_cluster = []
    for matrixConsensus in pMatrixNameConsensusList:
        difference_sample_consensus_list = []
        if pChromosomes is not None and len(pChromosomes) == 1:
            hic_ma_consensus = hm.hiCMatrix(pMatrixFile=pMatrixNameConsensus + '::' + matrixConsensus, pChrnameList=pChromosomes)
        else:
            hic_ma_consensus = hm.hiCMatrix(pMatrixFile=pMatrixNameConsensus + '::' + matrixConsensus)
            if pChromosomes:
                hic_ma_consensus.keepOnlyTheseChr(pChromosomes)

        for i, matrix in enumerate(pMatricesList):
            if pChromosomes is not None and len(pChromosomes) == 1:
                hic_ma_sample = hm.hiCMatrix(pMatrixFile=pMatrixName+ '::' +matrix, pChrnameList=pChromosomes)
            else:
                hic_ma_sample = hm.hiCMatrix(pMatrixFile=pMatrixName+ '::' +matrix)
                if pChromosomes:
                    hic_ma_sample.keepOnlyTheseChr(pChromosomes)
            difference_sample_consensus = (hic_ma_sample.matrix - hic_ma_consensus.matrix).sum()
            difference_sample_consensus_list.append(difference_sample_consensus)
        difference_sample_consensus_list_cluster.append(difference_sample_consensus_list)
    pQueue.put(difference_sample_consensus_list_cluster)
    return

def main(args=None):

    args = parse_arguments().parse_args(args)

    create_or_load_matrix = False

    # if args.createMatrix:
        
    # matrices_name_sample = 
    matrices_list_sample = cooler.fileops.list_coolers(args.matrix)
    matrices_list_consensus = cooler.fileops.list_coolers(args.consensusMatrix)
    threads = args.threads

    differences_per_cluster_threads =[None] * threads


    all_data_collected = False
    thread_done = [False] * threads
    log.debug('matrix read, starting processing')
    length_index = [None] * threads
    length_index[0] = 0
    matricesPerThread = len(matrices_list_sample) // threads
    queue = [None] * threads
    process = [None] * threads
    for i in range(threads):

        if i < threads - 1:
            matrices_name_list = matrices_list_sample[i * matricesPerThread:(i + 1) * matricesPerThread]
            # length_index[i + 1] = length_index[i] + len(matrices_name_list)
        else:
            matrices_name_list = matrices_list_sample[i * matricesPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=compute_consensus_sample_difference, kwargs=dict(
                            pMatrixName = args.matrix,
                            pMatricesList= matrices_name_list, 
                            pMatrixNameConsensus=args.consensusMatrix,
                            pMatrixNameConsensusList=matrices_list_consensus,
                            pChromosomes=args.chromosomes,
                            pQueue=queue[i]
            )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                differences_per_cluster_threads[i] = queue[i].get()
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

    log.debug('len(differences_per_cluster_threads) {}'.format(len(differences_per_cluster_threads)))
    log.debug('len(differences_per_cluster_threads[0]) {}'.format(len(differences_per_cluster_threads[0])))
    # log.debug('len(computed_results_threads) {}'.format(len(computed_results_threads)))

    differences_per_cluster = [None] * len(matrices_list_consensus) # --> number of clusters
    for computed_results_threads in differences_per_cluster_threads:
        for i, per_cluster in enumerate(computed_results_threads):
            if differences_per_cluster[i] is None:
                differences_per_cluster[i] = per_cluster
            else:
                differences_per_cluster[i] = np.concatenate((differences_per_cluster[i], per_cluster))
    log.debug('differences_per_cluster {}'.format(len(differences_per_cluster)))
    log.debug('differences_per_cluster[0] {}'.format(len(differences_per_cluster[0])))

    edge_count_matrix = csr_matrix((len(matrices_list_sample), len(matrices_list_sample)), dtype=int)
    edge_direction_matrix = csr_matrix((len(matrices_list_sample), len(matrices_list_sample)), dtype=int)

    for i, cluster in enumerate(differences_per_cluster):
        differences_per_cluster[i] = np.array(differences_per_cluster[i])
        argssorted_cluster = np.argsort(differences_per_cluster[i])
        j = 0
        while j < len(argssorted_cluster) - 1:
            edge_count_matrix[argssorted_cluster[j], argssorted_cluster[j+1]] += 1
            edge_direction_matrix[argssorted_cluster[j], argssorted_cluster[j+1]] += np.sign(argssorted_cluster[j])
            j += 1
    log.debug('computed edges')
    instances_edge_count, features_edge_count = edge_count_matrix.nonzero()
    # instances_edge_direction, features_edge_direction = edge_direction_matrix.nonzero()

    G=nx.DiGraph()
    for i in range(len(instances_edge_count)):
        if edge_direction_matrix.data[i] < 0:
            G.add_edge(instances_edge_count[i], features_edge_count[i], weight=edge_count_matrix[i])
        else:
            G.add_edge(features_edge_count[i], instances_edge_count[i], weight=edge_count_matrix[i])
    log.debug('computing hamiltonian path')
    # hamiltonian_path_output = nat.hamiltonian_path(G)

    # log.debug('hamiltonian_path {}'.format(hamiltonian_path_output))

    # nx.draw(G, node_size=10)
    # plt.savefig(args.outFileName, dpi=300)
    # plt.close()
    # nx.draw_random(G, node_size=10)
    # plt.savefig('random_'+args.outFileName, dpi=300)
    # plt.close()
    # nx.draw_circular(G, node_size=10)
    # plt.savefig('circular_'+args.outFileName, dpi=300)
    # plt.close()
    # nx.draw_spectral(G, node_size=10)
    # plt.savefig('spectral_'+args.outFileName, dpi=300)
    # plt.close()
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
