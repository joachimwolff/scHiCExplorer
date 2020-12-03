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

from holoviews.plotting.util import process_cmap

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
    parserOpt.add_argument('--dpi', '-d',
                           help='The dpi of the scatter plot.',
                           required=False,
                           default=300,
                           type=int)

    parserOpt.add_argument('--outFileName', '-o',
                           help='File name to save the resulting clusters',
                           required=True,
                           default='clusters.txt')
    parserOpt.add_argument('--cell_coloring_type', '-cct',
                           help='A two column list, first colum the cell names as stored in the scool file, second column the associated coloring for the scatter plot',
                           required=False)
    parserOpt.add_argument('--cell_coloring_batch', '-ccb',
                           help='A two column list, first colum the cell names as stored in the scool file, second column the associated coloring for the scatter plot',
                           required=False)
    parserOpt.add_argument('--latexTable', '-lt',
                           help='Return the overlap statistics if --cell_coloring_type is given as a latex table.')
    parserOpt.add_argument('--figuresize',
                           help='Fontsize in the plot for x and y axis.',
                           type=float,
                           nargs=2,
                           default=(15, 6),
                           metavar=('x-size', 'y-size'))
    parserOpt.add_argument('--colorMap',
                           help='Color map to use for the heatmap, supported are the categorical colormaps from holoviews: '
                           'http://holoviews.org/user_guide/Colormaps.html',
                           default='glasbey_dark')
    parserOpt.add_argument('--fontsize',
                           help='Fontsize in the plot for x and y axis.',
                           type=float,
                           default=15)
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module.',
                           required=False,
                           default=8,
                           type=int)
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def main(args=None):

    args = parse_arguments().parse_args(args)

    outputFolder = os.path.dirname(os.path.abspath(args.outFileName)) + '/'
    log.debug('outputFolder {}'.format(outputFolder))
    if args.cell_coloring_type:
        cell_name_cell_type_dict = {}

        cell_type_color_dict = {}
        color_cell_type_dict = {}
        cell_type_counter = 0
        with open(args.cell_coloring_type, 'r') as file:
            for i, line in enumerate(file.readlines()):
                line = line.strip()
                try:
                    cell_name, cell_type = line.split('\t')
                except Exception:
                    cell_name, cell_type = line.split('    ')
                cell_name_cell_type_dict[cell_name] = cell_type
                if cell_type not in cell_type_color_dict:
                    cell_type_color_dict[cell_type] = cell_type_counter
                    color_cell_type_dict[cell_type_counter] = cell_type
                    cell_type_counter += 1

    if args.cell_coloring_batch:
        cell_name_cell_type_dict_batch = {}

        cell_type_color_dict_batch = {}
        color_cell_type_dict_batch = {}
        cell_type_counter_batch = 0
        with open(args.cell_coloring_batch, 'r') as file:
            for i, line in enumerate(file.readlines()):
                line = line.strip()
                try:
                    cell_name, cell_type = line.split('\t')
                except Exception:
                    cell_name, cell_type = line.split('    ')
                cell_name_cell_type_dict_batch[cell_name] = cell_type
                if cell_type not in cell_type_color_dict_batch:
                    cell_type_color_dict_batch[cell_type] = cell_type_counter_batch
                    color_cell_type_dict_batch[cell_type_counter_batch] = cell_type
                    cell_type_counter_batch += 1

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
                                                       n_neighbors=reduce_to_dimension, affinity='nearest_neighbors', random_state=0, eigen_solver="arpack")

        labels_clustering = spectralClustering_object.fit_predict(neighborhood_matrix)
    elif args.clusterMethod == 'kmeans':
        kmeans_object = KMeans(n_clusters=args.numberOfClusters, random_state=0, n_jobs=args.threads, precompute_distances=True)
        labels_clustering = kmeans_object.fit_predict(neighborhood_matrix)

    if args.colorMap:
        colors = process_cmap(args.colorMap)
    if args.cell_coloring_type:
        if len(colors) < len(cell_type_color_dict):
            log.error('The chosen colormap offers too less values for the number of clusters.')
            exit(1)
        labels_clustering_cell_type = []
        for cell_name in matrices_list:
            labels_clustering_cell_type.append(cell_type_color_dict[cell_name_cell_type_dict[cell_name]])

        labels_clustering_cell_type = np.array(labels_clustering_cell_type)

        log.debug('labels_clustering_cell_type: {}'.format(len(labels_clustering_cell_type)))
        log.debug('matrices_list: {}'.format(len(matrices_list)))
    label_x = 'PC1'
    label_y = 'PC2'
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

        if args.cell_coloring_type:
            plt.figure(figsize=(args.figuresize[0], args.figuresize[1]))
            for i, color in enumerate(colors[:len(cell_type_color_dict)]):
                mask = labels_clustering_cell_type == i
                log.debug('plot cluster: {} {}'.format(color_cell_type_dict[i], np.sum(mask)))
                plt.scatter(neighborhood_matrix_knn[:, 0].T[mask], neighborhood_matrix_knn[:, 1].T[mask], color=color, label=str(color_cell_type_dict[i]), s=20, alpha=0.7)

            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=args.fontsize)
            plt.xticks([])
            plt.yticks([])
            plt.xlabel(label_x, fontsize=args.fontsize)
            plt.ylabel(label_y, fontsize=args.fontsize)
            if '.' not in args.createScatterPlot:
                args.createScatterPlot += '.png'
            scatter_plot_name = '.'.join(args.createScatterPlot.split('.')[:-1]) + '_cell_color.' + args.createScatterPlot.split('.')[-1]
            plt.tight_layout()
            plt.savefig(scatter_plot_name, dpi=args.dpi)
            plt.close()
        if args.cell_coloring_batch:
            if len(colors) < len(cell_type_color_dict_batch):
                log.error('The chosen colormap offers too less values for the number of clusters.')
                exit(1)
            labels_clustering_cell_type_batch = []
            for cell_name in matrices_list:
                labels_clustering_cell_type_batch.append(cell_type_color_dict_batch[cell_name_cell_type_dict_batch[cell_name]])

            labels_clustering_cell_type_batch = np.array(labels_clustering_cell_type_batch)

            log.debug('labels_clustering_cell_type: {}'.format(len(labels_clustering_cell_type_batch)))
            log.debug('matrices_list: {}'.format(len(matrices_list)))

            plt.figure(figsize=(args.figuresize[0], args.figuresize[1]))
            for i, color in enumerate(colors[:len(cell_type_color_dict_batch)]):
                mask = labels_clustering_cell_type_batch == i
                log.debug('plot cluster: {} {}'.format(color_cell_type_dict_batch[i], np.sum(mask)))
                plt.scatter(neighborhood_matrix_knn[:, 0].T[mask], neighborhood_matrix_knn[:, 1].T[mask], color=color, label=str(color_cell_type_dict_batch[i]), s=20, alpha=0.7)

            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=args.fontsize)
            plt.xticks([])
            plt.yticks([])
            plt.xlabel(label_x, fontsize=args.fontsize)
            plt.ylabel(label_y, fontsize=args.fontsize)
            if '.' not in args.createScatterPlot:
                args.createScatterPlot += '.png'
            scatter_plot_name = '.'.join(args.createScatterPlot.split('.')[:-1]) + '_cell_color_batch.' + args.createScatterPlot.split('.')[-1]
            plt.tight_layout()
            plt.savefig(scatter_plot_name, dpi=args.dpi)
            plt.close()

        plt.figure(figsize=(args.figuresize[0], args.figuresize[1]))
        for i, color in enumerate(colors[:args.numberOfClusters]):
            mask = labels_clustering == i
            plt.scatter(neighborhood_matrix_knn[:, 0].T[mask], neighborhood_matrix_knn[:, 1].T[mask], color=color, label=str(i), s=20, alpha=0.7)
        plt.legend(fontsize=args.fontsize)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=args.fontsize)

        plt.xticks([])
        plt.yticks([])
        plt.xlabel(label_x, fontsize=args.fontsize)
        plt.ylabel(label_y, fontsize=args.fontsize)
        if '.' not in args.createScatterPlot:
            args.createScatterPlot += '.png'
        scatter_plot_name = '.'.join(args.createScatterPlot.split('.')[:-1]) + '.' + args.createScatterPlot.split('.')[-1]
        plt.tight_layout()
        plt.savefig(scatter_plot_name, dpi=args.dpi)
        plt.close()

    if args.latexTable and args.cell_coloring_type:
        # compute overlap of cell_type find found clusters
        computed_clusters = set(labels_clustering)
        cell_type_amounts_dict = {}
        # percentage_threshold = 0.8

        for threshold in [0.7, 0.8, 0.9]:
            cell_type_amounts_dict[threshold] = {}
        with open(args.latexTable, 'w') as matches_file:
            header = '\\begin{table}[!htb]\n\\footnotesize\n\\begin{tabular}{|l'
            body = '\\hline Cluster '
            for i in range(len(color_cell_type_dict)):
                mask_cell_type = labels_clustering_cell_type == i
                header += '|c'
                body += '& ' + str(color_cell_type_dict[i]) + ' (' + str(np.sum(mask_cell_type)) + ' cells)'
            header += '|}\n'
            body += '\\\\\n'
            # body = ''
            for i in computed_clusters:
                body += '\\hline Cluster ' + str(i)
                mask_computed_clusters = labels_clustering == i
                body += ' (' + str(np.sum(mask_computed_clusters)) + ' cells)'
                for j in range(len(cell_type_color_dict)):
                    mask_cell_type = labels_clustering_cell_type == j
                    mask = mask_computed_clusters & mask_cell_type
                    number_of_matches = np.sum(mask)
                    body += '& ' + str(number_of_matches)

                    if number_of_matches != 1:
                        body += ' cells / '
                    else:
                        body += ' cell / '

                    body += '{:.2f}'.format((number_of_matches / np.sum(mask_computed_clusters)) * 100) + ' \\% '
                    for threshold in [0.7, 0.8, 0.9]:

                        if number_of_matches / np.sum(mask_computed_clusters) >= threshold:
                            if color_cell_type_dict[j] in cell_type_amounts_dict[threshold]:
                                cell_type_amounts_dict[threshold][color_cell_type_dict[j]] += number_of_matches
                            else:
                                cell_type_amounts_dict[threshold][color_cell_type_dict[j]] = number_of_matches
                        else:
                            if color_cell_type_dict[j] in cell_type_amounts_dict[threshold]:
                                continue
                            else:
                                cell_type_amounts_dict[threshold][color_cell_type_dict[j]] = 0
                body += '\\\\\n'
            body += '\\hline ' + '&' * len(cell_type_color_dict) + '\\\\\n'

            for threshold in [0.7, 0.8, 0.9]:
                body += '\\hline Correct identified $>{}\\%$'.format(int(threshold * 100))
                for i in range(len(cell_type_color_dict)):
                    mask_cell_type = labels_clustering_cell_type == i

                    if color_cell_type_dict[i] in cell_type_amounts_dict[threshold]:
                        body += '& ' + str(cell_type_amounts_dict[threshold][color_cell_type_dict[i]]) + ' / ' + str(np.sum(mask_cell_type)) + ' ('
                        body += '{:.2f}'.format((cell_type_amounts_dict[threshold][color_cell_type_dict[i]] / np.sum(mask_cell_type)) * 100)
                    else:
                        body += '& ' + str(0) + ' / ' + str(np.sum(mask_cell_type)) + ' ('
                        body += '{:.2f}'.format(0 / np.sum(mask_cell_type))

                    body += ' \\%)'
                body += '\\\\\n'
            body += '\\hline \n'
            body += '\\end{tabular}\n\\caption{}\n\\end{table}'

            matches_file.write(header)
            matches_file.write(body)

    matrices_cluster = list(zip(matrices_list, labels_clustering))
    np.savetxt(args.outFileName, matrices_cluster, fmt="%s")
