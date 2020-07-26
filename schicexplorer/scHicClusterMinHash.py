import argparse
import os
import time
from multiprocessing import Process, Queue


import numpy as np
from scipy.sparse import csr_matrix, vstack, save_npz, lil_matrix, dok_matrix
from sklearn.decomposition import PCA

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
    parserOpt.add_argument('--additionalPCA', '-pca',
                           help='Computes PCA on top of a k-nn. Can improve the cluster result.',
                           action='store_true')
    parserOpt.add_argument('--dimensionsPCA', '-dim_pca',
                           help='The number of dimensions from the PCA matrix that should be considered for clustering. Can improve the cluster result.',
                           default=20,
                           type=int)
    parserOpt.add_argument('--perChromosome', '-pc',
                           help='Computes the knn per chromosome and merge the different knns for clustering.',
                           action='store_true')
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

    raw_file_name = os.path.splitext(os.path.basename(args.outFileName))[0]
    if args.perChromosome:
        matrices_list = cell_name_list(args.matrix)

        neighborhood_matrix_knn = None
        cooler_obj = cooler.Cooler(args.matrix + '::' + matrices_list[0])
        # cooler_obj.chromnames

        if args.chromosomes is None:
            args.chromosomes = cooler_obj.chromnames
            log.debug('args.chromosomes was none')
        for chromosome in args.chromosomes:
            neighborhood_matrix, matrices_list = create_csr_matrix_all_cells(args.matrix, args.threads, [chromosome], outputFolder, raw_file_name, args.intraChromosomalContactsOnly)

            if len(neighborhood_matrix.data) == 0:
                log.debug('empty matrix chromosome {}'.format(chromosome))
                continue

            minHash_object = MinHash(n_neighbors=args.numberOfNearestNeighbors, number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
                                        shingle_size=4, fast=args.exactModeMinHash, maxFeatures=int(max(neighborhood_matrix.getnnz(1))), absolute_numbers=False)
            
            if args.shareOfMatrixToBeTransferred is not None and args.shareOfMatrixToBeTransferred > 0:
                if args.shareOfMatrixToBeTransferred > 1:
                    args.shareOfMatrixToBeTransferred = 1
                number_of_elements = neighborhood_matrix.shape[0]
                batch_size = int(np.floor(number_of_elements * args.shareOfMatrixToBeTransferred))
                if batch_size < 1:
                    batch_size = 1
                sub_matrix = neighborhood_matrix[0:batch_size, :]
                log.debug('chr {} len sub.data: {}'.format(chromosome, len(sub_matrix.data)))
                minHash_object.fit(neighborhood_matrix[0:batch_size, :])
                if batch_size < number_of_elements:
                    log.debug('partial fit')
                    for i in range(batch_size, neighborhood_matrix.shape[0], batch_size):
                        sub_matrix = neighborhood_matrix[i:i+batch_size, :]
                        log.debug('chr {} len sub.data: {}'.format(chromosome, len(sub_matrix.data)))
                        minHash_object.partial_fit(neighborhood_matrix[i:i+batch_size, :])
            else:
                minHash_object.fit(neighborhood_matrix)
            log.debug('147')
            if neighborhood_matrix_knn is None:
                neighborhood_matrix_knn = minHash_object.kneighbors_graph(mode='distance')
            else:
                neighborhood_matrix_knn += minHash_object.kneighbors_graph(mode='distance')
            del minHash_object
        if args.dimensionsPCA:
            pca = PCA(n_components = min(neighborhood_matrix_knn.shape) - 1)
            neighborhood_matrix_knn = pca.fit_transform(neighborhood_matrix_knn.todense())
            if args.dimensionsPCA:
                args.dimensionsPCA = min(args.dimensionsPCA, neighborhood_matrix_knn.shape[0])
                neighborhood_matrix_knn = neighborhood_matrix_knn[:, :args.dimensionsPCA]
        
        if args.clusterMethod == 'spectral':
            spectralClustering_object = SpectralClustering(n_clusters=args.numberOfClusters, n_jobs=args.threads,
                                                        n_neighbors=reduce_to_dimension, affinity='nearest_neighbors', random_state=0)

            labels_clustering = spectralClustering_object.fit_predict(neighborhood_matrix_knn)
        elif args.clusterMethod == 'kmeans':
            kmeans_object = KMeans(n_clusters=args.numberOfClusters, random_state=0, n_jobs=args.threads, precompute_distances=True)
            labels_clustering = kmeans_object.fit_predict(neighborhood_matrix_knn)


        # minHash_object.fit()
    else:
        
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

            # if args.pca:

            labels_clustering = minHashClustering.fit_predict(neighborhood_matrix, pSaveMemory=args.shareOfMatrixToBeTransferred, pPca=args.additionalPCA, pPcaDimensions=args.dimensionsPCA)
            log.debug('create label matrix assoziation')
        elif args.clusterMethod == 'kmeans':
            log.debug('kmeans clustering')
            kmeans_object = KMeans(n_clusters=args.numberOfClusters, random_state=0, n_jobs=args.threads, precompute_distances=True)
            minHash_object = MinHash(n_neighbors=args.numberOfNearestNeighbors, number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
                                    shingle_size=4, fast=args.exactModeMinHash, maxFeatures=int(max(neighborhood_matrix.getnnz(1))), absolute_numbers=False)
            minHashClustering = MinHashClustering(minHashObject=minHash_object, clusteringObject=kmeans_object)
            log.debug('kmeans clustering fit predict')

            labels_clustering = minHashClustering.fit_predict(neighborhood_matrix, pSaveMemory=args.shareOfMatrixToBeTransferred, pPca=args.additionalPCA, pPcaDimensions=args.dimensionsPCA)


    matrices_cluster = list(zip(matrices_list, labels_clustering))
    np.savetxt(args.outFileName, matrices_cluster, fmt="%s")
