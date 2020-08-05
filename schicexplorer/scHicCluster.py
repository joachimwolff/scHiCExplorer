import argparse
import os
from multiprocessing import Process, Queue
import time


import logging
log = logging.getLogger(__name__)
from scipy import linalg
import cooler
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
from hicmatrix import HiCMatrix as hm

import numpy as np
from scipy.sparse import csr_matrix

from schicexplorer._version import __version__
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
    parserOpt.add_argument('--intraChromosomalContactsOnly', '-ic',
                           help='This option loads only the intra-chromosomal contacts. Can improve the cluster result if data is very noisy.',
                           action='store_true')
    parserOpt.add_argument('--additionalPCA', '-pca',
                           help='Computes PCA on top of a k-nn. Can improve the cluster result.',
                           action='store_true')
    parserOpt.add_argument('--dimensionsPCA', '-dim_pca',
                           help='The number of dimensions from the PCA matrix that should be considered for clustering. Can improve the cluster result.',
                           default=20,
                           type=int)
    parserOpt.add_argument('--dimensionReductionMethod', '-drm',
                           help='Dimension reduction methods, knn with euclidean distance, pca',
                           choices=['none', 'knn', 'pca'],
                           default='none')
    parserOpt.add_argument('--createScatterPlot', '-csp',
                           help='Create a scatter plot for the clustering, the x and y are the first and second principal component of the computed k-nn graph.',
                           required=False,
                           default=None)
    parserOpt.add_argument('--numberOfNearestNeighbors', '-k',
                           help='Number of to be used computed nearest neighbors for the knn graph. Default is either the default value or the number of the provided cells, whatever is smaller.',
                           required=False,
                           default=100,
                           type=int)
    parserOpt.add_argument('--colorMap',
                           help='Color map to use for the heatmap. Available '
                           'values can be seen here: '
                           'http://matplotlib.org/examples/color/colormaps_reference.html',
                           default='tab20')
    parserOpt.add_argument('--dpi', '-d',
                           help='The dpi of the scatter plot.',
                           required=False,
                           default=300,
                           type=int)
    parserOpt.add_argument('--fontsize',
                           help='Fontsize in the plot for x and y axis.',
                           type=float,
                           default=10)
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

    outputFolder = os.path.dirname(os.path.abspath(args.outFileName)) + '/'
    log.debug('outputFolder {}'.format(outputFolder))

    raw_file_name = os.path.splitext(os.path.basename(args.outFileName))[0]
    neighborhood_matrix, matrices_list = create_csr_matrix_all_cells(args.matrix, args.threads, args.chromosomes, outputFolder, raw_file_name, args.intraChromosomalContactsOnly)

    reduce_to_dimension = neighborhood_matrix.shape[0] - 1
    if args.dimensionReductionMethod == 'knn':

        if args.numberOfNearestNeighbors > reduce_to_dimension:
            args.numberOfNearestNeighbors = reduce_to_dimension
        nbrs = NearestNeighbors(n_neighbors=args.numberOfNearestNeighbors, algorithm='ball_tree', n_jobs=args.threads).fit(neighborhood_matrix)
        neighborhood_matrix = nbrs.kneighbors_graph(mode='distance')

        if args.additionalPCA:
            pca = PCA(n_components=min(neighborhood_matrix.shape) - 1)
            neighborhood_matrix = pca.fit_transform(neighborhood_matrix.todense())
            if args.dimensionsPCA:
                args.dimensionsPCA = min(args.dimensionsPCA, neighborhood_matrix.shape[0])
                neighborhood_matrix = neighborhood_matrix[:, :args.dimensionsPCA]
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

    if args.createScatterPlot:
        if args.dimensionReductionMethod == 'none':
            log.warning('Raw matrix clustering scatter plot needs to compute a PCA and can request large amount (> 100 GB) of memory.')

        log.debug('args.additionalPCA  {}'.format(args.additionalPCA))
        log.debug('args.dimensionReductionMethod  {}'.format(args.dimensionReductionMethod))

        if args.dimensionReductionMethod == 'none' or (args.dimensionReductionMethod == 'knn' and not args.additionalPCA):
            log.debug('compute pca')

            pca = PCA(n_components=min(neighborhood_matrix.shape) - 1)
            neighborhood_matrix_knn = pca.fit_transform(neighborhood_matrix.todense())
            log.debug('compute pca')
        else:
            log.debug('already computed pca')

            neighborhood_matrix_knn = neighborhood_matrix
        plt.figure(figsize=(15, 8))

        list(set(labels_clustering))
        plt.figure(figsize=(15, 8), dpi=80)
        cmap = get_cmap(args.colorMap)
        colors = cmap.colors
        for i, color in enumerate(colors[:args.numberOfClusters]):
            mask = labels_clustering == i
            plt.scatter(neighborhood_matrix_knn[:, 0].T[mask], neighborhood_matrix_knn[:, 1].T[mask], color=color, label=str(i), s=50, alpha=0.7)
        plt.legend(fontsize=args.fontsize)
        plt.xticks([])
        plt.yticks([])
        plt.xlabel('PC1', fontsize=args.fontsize)
        plt.ylabel('PC2', fontsize=args.fontsize)
        log.debug('args.createScatterPlot {}'.format(args.createScatterPlot))
        if '.' not in args.createScatterPlot:
            args.createScatterPlot += '.png'
        scatter_plot_name = '.'.join(args.createScatterPlot.split('.')[:-1]) + '_pc1_pc2.' + args.createScatterPlot.split('.')[-1]
        log.debug('scatter_plot_name {}'.format(scatter_plot_name))
        plt.savefig(scatter_plot_name, dpi=args.dpi)
        plt.close()

        for i, color in enumerate(colors[:args.numberOfClusters]):
            mask = labels_clustering == i
            plt.scatter(neighborhood_matrix_knn[:, 1].T[mask], neighborhood_matrix_knn[:, 2].T[mask], color=color, label=str(i), s=50, alpha=0.7)
        plt.legend(fontsize=args.fontsize)
        plt.xticks([])
        plt.yticks([])
        plt.xlabel('PC2', fontsize=args.fontsize)
        plt.ylabel('PC3', fontsize=args.fontsize)
        scatter_plot_name = '.'.join(args.createScatterPlot.split('.')[:-1]) + '_pc2_pc3.' + args.createScatterPlot.split('.')[-1]
        plt.savefig(scatter_plot_name, dpi=args.dpi)
        plt.close()

    matrices_cluster = list(zip(matrices_list, labels_clustering))
    np.savetxt(args.outFileName, matrices_cluster, fmt="%s")
