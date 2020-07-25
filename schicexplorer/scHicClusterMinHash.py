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
    parserOpt.add_argument('--intraChromosomalContactsOnly', '-ic',
                           help='This option loads only the intra-chromosomal contacts. Can improve the cluster result if data is very noisy.',
                           action='store_true')
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



def main(args=None):

    args = parse_arguments().parse_args(args)

    outputFolder = os.path.dirname(os.path.abspath(args.outFileName)) + '/'
    log.debug('outputFolder {}'.format(outputFolder))

    raw_file_name = os.path.splitext(os.path.basename(args.outFileName))[0]
    neighborhood_matrix, matrices_list = create_csr_matrix_all_cells(args.matrix, args.threads, args.chromosomes, outputFolder, raw_file_name, args.intraChromosomalContactsOnly)
    if args.saveIntermediateRawMatrix:
        save_npz(args.saveIntermediateRawMatrix, neighborhood_matrix)
    if args.clusterMethod == 'spectral':
        log.debug('spectral clustering')
        spectral_object = SpectralClustering(n_clusters=args.numberOfClusters, affinity='nearest_neighbors', n_jobs=args.threads, random_state=0)
        log.debug('spectral clustering fit predict')
        minHash_object = MinHash(n_neighbors=args.numberOfNearestNeighbors, number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
                                 shingle_size=4, fast=args.exactModeMinHash, maxFeatures=int(max(neighborhood_matrix.getnnz(1))), absolute_numbers=False)
        minHashClustering = MinHashClustering(minHashObject=minHash_object, clusteringObject=spectral_object)

        labels_clustering = minHashClustering.fit_predict(neighborhood_matrix, pSaveMemory=args.shareOfMatrixToBeTransferred)
        log.debug('create label matrix assoziation')
    elif args.clusterMethod == 'kmeans':
        log.debug('kmeans clustering')
        kmeans_object = KMeans(n_clusters=args.numberOfClusters, random_state=0, n_jobs=args.threads, precompute_distances=True)
        minHash_object = MinHash(n_neighbors=args.numberOfNearestNeighbors, number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
                                 shingle_size=4, fast=args.exactModeMinHash, maxFeatures=int(max(neighborhood_matrix.getnnz(1))), absolute_numbers=False)
        minHashClustering = MinHashClustering(minHashObject=minHash_object, clusteringObject=kmeans_object)
        log.debug('kmeans clustering fit predict')

        labels_clustering = minHashClustering.fit_predict(neighborhood_matrix, pSaveMemory=args.shareOfMatrixToBeTransferred)

    matrices_cluster = list(zip(matrices_list, labels_clustering))
    np.savetxt(args.outFileName, matrices_cluster, fmt="%s")
