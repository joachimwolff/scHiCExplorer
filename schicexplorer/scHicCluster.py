import argparse
import os
from multiprocessing import Process, Queue
import time


import logging
log = logging.getLogger(__name__)
from scipy import linalg
import cooler

from sklearn.cluster import KMeans, SpectralClustering
from sklearn.neighbors import NearestNeighbors
from hicmatrix import HiCMatrix as hm

import numpy as np
from scipy.sparse import csr_matrix

from schicexplorer._version import __version__
# from schicexplorer.utilities import cell_name_list
from schicexplorer.utilities import cell_name_list, create_csr_matrix_all_cells


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description='scHicCluster uses kmeans or spectral clustering to associate each cell to a cluster and therefore to its cell cycle. '
        'The clustering can be run on the raw data, on a kNN computed via the exact euclidean distance or via PCA. '
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
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--chromosomes',
                           help='List of to be plotted chromosomes',
                           nargs='+')

    parserOpt.add_argument('--dimensionReductionMethod', '-drm',
                           help='Dimension reduction methods, knn with euclidean distance, pca',
                           choices=['none', 'knn', 'pca'],
                           default='none')
    parserOpt.add_argument('--numberOfNearestNeighbors', '-k',
                           help='Number of to be used computed nearest neighbors for the knn graph. Default is either the default value or the number of the provided cells, whatever is smaller.',
                           required=False,
                           default=100,
                           type=int)
    parserOpt.add_argument('--outFileName', '-o',
                           help='File name to save the resulting clusters',
                           required=True,
                           default='clusters.txt')
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

  
    neighborhood_matrix, matrices_list = create_csr_matrix_all_cells(args.matrix, args.threads, args.chromosomes)

    reduce_to_dimension = neighborhood_matrix.shape[0] - 1
    if args.dimensionReductionMethod == 'knn':

        if args.numberOfNearestNeighbors > reduce_to_dimension:
            args.numberOfNearestNeighbors = reduce_to_dimension
        nbrs = NearestNeighbors(n_neighbors=args.numberOfNearestNeighbors, algorithm='ball_tree', n_jobs=args.threads).fit(neighborhood_matrix)
        neighborhood_matrix = nbrs.kneighbors_graph(mode='distance')

    elif args.dimensionReductionMethod == 'pca':
        corrmatrix = np.cov(neighborhood_matrix.todense())
        evals, eigs = linalg.eig(corrmatrix)
        neighborhood_matrix = eigs[:, :reduce_to_dimension].transpose()

    if args.clusterMethod == 'spectral':
        spectralClustering_object = SpectralClustering(n_clusters=args.numberOfClusters, n_jobs=args.threads,
                                                       n_neighbors=reduce_to_dimension, affinity='nearest_neighbors', random_state=0)

        labels_clustering = spectralClustering_object.fit_predict(neighborhood_matrix)
    elif args.clusterMethod == 'kmeans':
        kmeans_object = KMeans(n_clusters=args.numberOfClusters, random_state=0, n_jobs=args.threads, precompute_distances=True)
        labels_clustering = kmeans_object.fit_predict(neighborhood_matrix)

    matrices_cluster = list(zip(matrices_list, labels_clustering))
    np.savetxt(args.outFileName, matrices_cluster, fmt="%s")
