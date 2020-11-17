import argparse
import os
import time
from multiprocessing import Process, Queue


import numpy as np
from scipy.sparse import csr_matrix, vstack, save_npz, lil_matrix, dok_matrix, isspmatrix_csr
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import logging
log = logging.getLogger(__name__)

import cooler
from sklearn.cluster import SpectralClustering, KMeans, AgglomerativeClustering, Birch

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from sparse_neighbors_search import MinHash
from sparse_neighbors_search import MinHashClustering

from hicmatrix import HiCMatrix as hm

from schicexplorer._version import __version__
from schicexplorer.utilities import cell_name_list, create_csr_matrix_all_cells, open_and_store_matrix

from holoviews.plotting.util import process_cmap

def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description='scHicClusterMinHash uses kmeans, spectral, agglomerative (ward, complete, average or single) or birch clustering to associate each cell to a cluster and therefore to its cell cycle. '
        'The clustering is applied on dimension reduced data based on an approximate kNN search with the local sensitive hashing technique MinHash. This approach reduces the number of dimensions from samples * (number of bins)^2 to samples * samples. '
        'Please consider also the other clustering and dimension reduction approaches of the scHicExplorer suite. They can give you better results, '
        'can be faster or less memory demanding.'

    )

    parserRequired = parser.add_argument_group('Required arguments')

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
                                choices=['spectral', 'kmeans', 'agglomerative_ward', 'agglomerative_complete', 'agglomerative_average', 'agglomerative_single', 'birch'],
                                default='spectral')
    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the resulting clusters',
                                required=True,
                                default='clusters.txt')
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--euclideanModeMinHash', '-em',
                           help='This option uses the number of hash collisions is only for a candidate set selection and computes on them the euclidean distance.',
                           action='store_false')
    parserOpt.add_argument('--intraChromosomalContactsOnly', '-ic',
                           help='This option loads only the intra-chromosomal contacts. Can improve the cluster result if data is very noisy.',
                           action='store_true')
    parserOpt.add_argument('--saveIntermediateRawMatrix', '-sirm',
                           help='This option activates the save of the intermediate raw scHi-C matrix.',
                           required=False)
    parserOpt.add_argument('--createScatterPlot', '-csp',
                           help='Create a scatter plot for the clustering, the x and y are the first and second principal component of the computed k-nn graph.',
                           required=False,
                           default='scatterPlot.eps')
    parserOpt.add_argument('--dpi',
                           help='The dpi of the scatter plot.',
                           required=False,
                           default=300,
                           type=int)
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
                           help='Which share of rows shall be transferred from Python to C++ at once. If `--saveMemory` is active, value is interpreted as the share to loaded at once to memory. Values between 0 and 1, the more are transferred at once, the larger the memory usage is. The less rows are transferred, the slower the computation is.',
                           required=False,
                           default=0.25,
                           type=float)
    parserOpt.add_argument('--saveMemory', '-sm',
                           help='Load data only with one core, this method saves memory but is significantly slower.',
                           required=False,
                           action='store_true')
    parserOpt.add_argument('--cell_coloring_type', '-cct',
                            help='A two column list, first colum the cell names as stored in the scool file, second column the associated coloring for the scatter plot',
                            required=False)
    parserOpt.add_argument('--cell_coloring_batch', '-ccb',
                            help='A two column list, first colum the cell names as stored in the scool file, second column the associated coloring for the scatter plot',
                            required=False)
    parserOpt.add_argument('--noPCA',
                           help='Do not computes PCA on top of a k-nn.',
                           action='store_true')
    parserOpt.add_argument('--noUMAP',
                           help='Do not computes UMP on top of a k-nn/PCA.',
                           action='store_true')
    parserOpt.add_argument('--dimensionsPCA', '-dim_pca',
                           help='The number of dimensions from the PCA matrix that should be considered for clustering. Can improve the cluster result.',
                           default=100,
                           type=int)
    parserOpt.add_argument('--colorMap',
                           help='Color map to use for the heatmap. Available '
                           'values can be seen here: '
                           'http://matplotlib.org/examples/color/colormaps_reference.html',
                           default='tab10')
    parserOpt.add_argument('--fontsize',
                           help='Fontsize in the plot for x and y axis.',
                           type=float,
                           default=15)
    parserOpt.add_argument('--distance', '-d',
                           help='Contact distance to consider',
                           type=float,
                           default=None)
    parserOpt.add_argument('--figuresize',
                           help='Fontsize in the plot for x and y axis.',
                           type=float,
                           nargs=2,
                           default=(15, 6),
                           metavar=('x-size', 'y-size'))
    parserOpt.add_argument('--chromosomes',
                           help='List of to be computed chromosomes',
                           nargs='+')
    parserOpt.add_argument('--perChromosome', '-pc',
                           help='Computes the knn per chromosome and merge the different knns for clustering.',
                           action='store_true')
    parserOpt.add_argument('--absoluteValues', '-av',
                           help='Return the number of hash collisions as measure instead of 0 - 1 normalized values.',
                           action='store_true')
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module.',
                           required=False,
                           default=8,
                           type=int)
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    parserOptUmap = parser.add_argument_group('Optional umap arguments', 'Arguments for umap embedding. Please consider its documentation for details: https://umap-learn.readthedocs.io/en/latest/api.html#umap.umap_.UMAP')
    parserOptUmap.add_argument('--umap_n_neighbors',
                           help='Number of neighbors',
                           type=int,
                           default=30)
    parserOptUmap.add_argument('--umap_n_components',
                           help='Number of components',
                           type=int,
                           default=2)
    parserOptUmap.add_argument('--umap_metric',
                           help='Metric of umap.',
                           type=str,
                           default='canberra')
    parserOptUmap.add_argument('--umap_n_epochs',
                           help='Number of epochs',
                           type=int,
                           default=None)
    parserOptUmap.add_argument('--umap_learning_rate',
                           help='Learning rate',
                           type=float,
                           default=1.0)
    parserOptUmap.add_argument('--umap_init',
                           help='Initialization method',
                           type=str,
                           default='spectral')
    parserOptUmap.add_argument('--umap_min_dist',
                           help='Minimum distance of two neighbors',
                           type=float,
                           default=0.3)
    parserOptUmap.add_argument('--umap_spread',
                           help='Spread',
                           type=float,
                           default=1.0)
    parserOptUmap.add_argument('--umap_set_op_mix_ratio',
                           help='set_op_mix_ratio',
                           type=float,
                           default=1.0)
    parserOptUmap.add_argument('--umap_local_connectivity',
                           help='local connectivity',
                           type=float,
                           default=1.0)
    parserOptUmap.add_argument('--umap_repulsion_strength',
                           help='repulsion strength',
                           type=float,
                           default=1.0)
    parserOptUmap.add_argument('--umap_negative_sample_rate',
                           help='negative sample rate',
                           type=int,
                           default=5)
    parserOptUmap.add_argument('--umap_transform_queue_size',
                           help='transform queue size',
                           type=float,
                           default=4.0)
    parserOptUmap.add_argument('--umap_a',
                           help='a',
                           type=float,
                           default=None)
    parserOptUmap.add_argument('--umap_b',
                           help='b',
                           type=float,
                           default=None)
    parserOptUmap.add_argument('--umap_angular_rp_forest',
                           help='angular rp forest',
                           action='store_true')
    parserOptUmap.add_argument('--umap_target_n_neighbors',
                           help='target number of neighbors',
                           type=int,
                           default=-1)
    parserOptUmap.add_argument('--umap_target_metric',
                           help='target metric',
                           type=str,
                           default='categorical')
    parserOptUmap.add_argument('--umap_target_weight',
                           help='target weight',
                           type=float,
                           default=0.5)
    parserOptUmap.add_argument('--umap_force_approximation_algorithm',
                           help='force approximation algorithm',
                           action='store_true')
    parserOptUmap.add_argument('--umap_verbose',
                           help='verbose',
                           action='store_true')
    parserOptUmap.add_argument('--umap_unique',
                           help='Contact distance to consider',
                           action='store_true')

    return parser


def main(args=None):

    args = parse_arguments().parse_args(args)
    # if args.threads <= 4:
    #     log.error('')
    #     exit(1)
    outputFolder = os.path.dirname(os.path.abspath(args.outFileName)) + '/'

    raw_file_name = os.path.splitext(os.path.basename(args.outFileName))[0]

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
                except:
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
                except:
                    cell_name, cell_type = line.split('    ')
                cell_name_cell_type_dict_batch[cell_name] = cell_type
                if cell_type not in cell_type_color_dict_batch:
                    cell_type_color_dict_batch[cell_type] = cell_type_counter_batch
                    color_cell_type_dict_batch[cell_type_counter_batch] = cell_type
                    cell_type_counter_batch += 1
        # log.debug('cell_name_cell_type_dict,keys() {}'.format(cell_name_cell_type_dict.keys()))
        # log.debug('len(cell_type_color_dict) {}'.format(len(cell_type_color_dict)))
    if args.clusterMethod == 'spectral':
        cluster_object = SpectralClustering(n_clusters=args.numberOfClusters, affinity='nearest_neighbors', n_jobs=args.threads, random_state=0)
    elif args.clusterMethod == 'kmeans':
        cluster_object = KMeans(n_clusters=args.numberOfClusters, random_state=0, n_jobs=args.threads, precompute_distances=True)
    elif args.clusterMethod.startswith('agglomerative'):
        for linkage in ['ward', 'complete', 'average', 'single']:
            if linkage in args.clusterMethod:
                cluster_object = AgglomerativeClustering(n_clusters=args.numberOfClusters, linkage=linkage)
                break
    elif args.clusterMethod == 'birch':
        cluster_object = Birch(n_clusters=args.numberOfClusters)
    else:
        log.error('No valid cluster method given: {}'.format(args.clusterMethod))

    umap_params_dict = {}

    if not args.noUMAP:
        for param in vars(args):
            if 'umap_' in param:
                 umap_params_dict[param] = vars(args)[param]
    # log.debug(umap_params_dict)

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

            # minHash_object = MinHash(n_neighbors=args.numberOfNearestNeighbors, number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
            #                             shingle_size=4, fast=args.exactModeMinHash, maxFeatures=int(max(neighborhood_matrix.getnnz(1))), absolute_numbers=False)

            minHash_object = MinHash(n_neighbors=args.numberOfNearestNeighbors, number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
                                    shingle_size=5, fast=args.euclideanModeMinHash, maxFeatures=int(max(neighborhood_matrix.getnnz(1))), absolute_numbers=False, max_bin_size=100000,
                                    minimal_blocks_in_common=400, excess_factor=1, prune_inverse_index=False)
            
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

            # compute knn per chromosome
            neighborhood_matrix_knn_pca = minHash_object.kneighbors_graph(mode='distance')
            if args.dimensionsPCA:
                pca = PCA(n_components = min(neighborhood_matrix_knn_pca.shape) - 1)
                neighborhood_matrix_knn_pca = pca.fit_transform(neighborhood_matrix_knn_pca.todense())
                if args.dimensionsPCA:
                    args.dimensionsPCA = min(args.dimensionsPCA, neighborhood_matrix_knn_pca.shape[0])
                    neighborhood_matrix_knn_pca = neighborhood_matrix_knn_pca[:, :args.dimensionsPCA]
                
                    # log.debug('neighborhood_matrix_knn_pca {}'.format(neighborhood_matrix_knn_pca))
                    log.debug('neighborhood_matrix_knn_pca shape{}'.format(neighborhood_matrix_knn_pca.shape))

            if neighborhood_matrix_knn is None:
                neighborhood_matrix_knn = neighborhood_matrix_knn_pca#minHash_object.kneighbors_graph(mode='distance')
            else:
                neighborhood_matrix_knn += neighborhood_matrix_knn_pca
                # neighborhood_matrix_knn = np.hstack((neighborhood_matrix_knn, neighborhood_matrix_knn_pca))#minHash_object.kneighbors_graph(mode='distance')
            log.debug('neighborhood_matrix_knn shape{}'.format(neighborhood_matrix_knn.shape))
            del minHash_object

        # if args.dimensionsPCA:
        #     pca = PCA(n_components = min(neighborhood_matrix_knn.shape) - 1)
        #     neighborhood_matrix_knn = pca.fit_transform(neighborhood_matrix_knn)
        #     if args.dimensionsPCA:
        #         log.debug('neighborhood_matrix_knn All before pca shape{}'.format(neighborhood_matrix_knn.shape))

        #         args.dimensionsPCA = min(args.dimensionsPCA, neighborhood_matrix_knn.shape[0])
        #         neighborhood_matrix_knn = neighborhood_matrix_knn[:, :args.dimensionsPCA]
        #         log.debug('neighborhood_matrix_knn All after pca shape{}'.format(neighborhood_matrix_knn.shape))
        if args.clusterMethod == 'spectral':
            spectralClustering_object = SpectralClustering(n_clusters=args.numberOfClusters, n_jobs=args.threads,
                                                        n_neighbors=reduce_to_dimension, affinity='nearest_neighbors', random_state=0)

            labels_clustering = spectralClustering_object.fit_predict(neighborhood_matrix_knn)
        elif args.clusterMethod == 'kmeans':
            kmeans_object = KMeans(n_clusters=args.numberOfClusters, random_state=0, n_jobs=args.threads, precompute_distances=True)
            labels_clustering = kmeans_object.fit_predict(neighborhood_matrix_knn)
        # minHashClustering = MinHashClustering(minHashObject=None, clusteringObject=cluster_object)
        # minHashClustering._precomputed_graph = neighborhood_matrix_knn
        # labels_clustering = minHashClustering.predict(minHashClustering._precomputed_graph, pPca=args.noPCA, pPcaDimensions=args.dimensionsPCA)


        # minHash_object.fit()
    else:
        if args.saveMemory:
            matrices_list = cell_name_list(args.matrix)
            max_nnz = 0
            for matrix in matrices_list:
                cooler_obj = cooler.Cooler(args.matrix + '::' + matrix)
                nnz = cooler_obj.info['nnz']
                if max_nnz < nnz:
                    max_nnz = nnz
            minHash_object = None
            matricesPerRun = int(len(matrices_list) * args.shareOfMatrixToBeTransferred)
            if matricesPerRun < 1:
                matricesPerRun = 1
            chromosome_indices = None
            if args.intraChromosomalContactsOnly:
                cooler_obj = cooler.Cooler(args.matrix + '::' + matrices_list[0])
                binsDataFrame = cooler_obj.bins()[:]
                chromosome_indices = {}
                for chromosome in cooler_obj.chromnames:
                    chromosome_indices[chromosome] = np.array(binsDataFrame.index[binsDataFrame['chrom'] == chromosome].tolist())

            for j, i in enumerate(range(0, len(matrices_list), matricesPerRun)):
                if i < len(matrices_list) - 1:
                    matrices_share = matrices_list[i:i + matricesPerRun]
                else:
                    matrices_share = matrices_list[i:]
                neighborhood_matrix, matrices_list_share = open_and_store_matrix(args.matrix, matrices_share, 0, len(matrices_share),
                                                                                args.chromosomes, args.intraChromosomalContactsOnly, chromosome_indices)
                if minHash_object is None:
                    minHash_object = MinHash(n_neighbors=args.numberOfNearestNeighbors, number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
                                            shingle_size=0, fast=args.euclideanModeMinHash, maxFeatures=int(max_nnz), absolute_numbers=False)

                if j == 0:
                    minHash_object.fit(neighborhood_matrix)
                else:
                    minHash_object.partial_fit(X=neighborhood_matrix)

            precomputed_graph = minHash_object.kneighbors_graph(mode='distance')

            if not args.noPCA:

                pca = PCA(n_components=min(precomputed_graph.shape) - 1)
                precomputed_graph = pca.fit_transform(precomputed_graph.todense())

                if args.dimensionsPCA:
                    args.dimensionsPCA = min(args.dimensionsPCA, precomputed_graph.shape[0])
                    cluster_object.fit(precomputed_graph[:, :args.dimensionsPCA])
                    # return
            else:
                try:
                    cluster_object.fit(precomputed_graph)
                except Exception:
                    cluster_object.fit(precomputed_graph.todense())
            minHashClustering = MinHashClustering(minHashObject=minHash_object, clusteringObject=cluster_object)
            minHashClustering._precomputed_graph = precomputed_graph

        else:
            neighborhood_matrix, matrices_list = create_csr_matrix_all_cells(args.matrix, args.threads, args.chromosomes, outputFolder, raw_file_name, args.intraChromosomalContactsOnly, pDistance=args.distance)

            if args.saveIntermediateRawMatrix:
                save_npz(args.saveIntermediateRawMatrix, neighborhood_matrix)

        if not args.saveMemory:
            minHash_object = MinHash(n_neighbors=args.numberOfNearestNeighbors, number_of_hash_functions=args.numberOfHashFunctions, number_of_cores=args.threads,
                                    shingle_size=5, fast=args.euclideanModeMinHash, maxFeatures=int(max(neighborhood_matrix.getnnz(1))), absolute_numbers=False, max_bin_size=100000,
                                    minimal_blocks_in_common=100, excess_factor=1, prune_inverse_index=False)
            minHashClustering = MinHashClustering(minHashObject=minHash_object, clusteringObject=cluster_object)
            minHashClustering.fit(X=neighborhood_matrix, pSaveMemory=args.shareOfMatrixToBeTransferred, pPca=(not args.noPCA), pPcaDimensions=args.dimensionsPCA, pUmap=(not args.noUMAP), pUmapDict=umap_params_dict)
            # log.debug('251 type(){}'.format(type(minHashClustering._precomputed_graph)))

        if args.noPCA and args.noUMAP:
            mask = np.isnan(minHashClustering._precomputed_graph.data)
            minHashClustering._precomputed_graph.data[mask] = 0

            mask = np.isinf(minHashClustering._precomputed_graph.data)
            minHashClustering._precomputed_graph.data[mask] = 0

        labels_clustering = minHashClustering.predict(minHashClustering._precomputed_graph, pPca=args.noPCA, pPcaDimensions=args.dimensionsPCA)


    if args.createScatterPlot:
        if args.noPCA  and args.noUMAP:
            pca = PCA(n_components=min(minHashClustering._precomputed_graph.shape) - 1)
            neighborhood_matrix_knn = pca.fit_transform(minHashClustering._precomputed_graph.todense())
        else:
            if not args.perChromosome:
                neighborhood_matrix_knn = minHashClustering._precomputed_graph
        # neighborhood_matrix_knn = minHashClustering._precomputed_graph
       

        list(set(labels_clustering))
        # cmap = get_cmap(args.colorMap)
        # colors = cmap.colors

        colors = process_cmap(args.colorMap)
        
        try:
            neighborhood_matrix_knn = neighborhood_matrix_knn.toarray()
        except Exception:
            pass
        
        # if not (args.noPCA):
        label_x = 'PC1'
        label_y = 'PC2'
        if not (args.noUMAP):
            label_x = 'UMAP1'
            label_y = 'UMAP2'
        if args.cell_coloring_type:
            if len(colors) < len(cell_type_color_dict):
                log.error('The chosen colormap offers too less values for the number of clusters.')
                exit(1)
            labels_clustering_cell_type = []
            for cell_name in matrices_list:
                # if cell_name in cell_name_cell_type_dict:
                #     if cell_name_cell_type_dict[cell_name] in cell_type_color_dict:
                labels_clustering_cell_type.append(cell_type_color_dict[cell_name_cell_type_dict[cell_name]])

            labels_clustering_cell_type = np.array(labels_clustering_cell_type)

            log.debug('labels_clustering_cell_type: {}'.format(len(labels_clustering_cell_type)))
            log.debug('matrices_list: {}'.format(len(matrices_list)))


            plt.figure(figsize=(args.figuresize[0], args.figuresize[1]))
            for i, color in enumerate(colors[:len(cell_type_color_dict)]):
                mask = labels_clustering_cell_type == i
                log.debug('plot cluster: {} {}'.format(color_cell_type_dict[i], np.sum(mask)))
                plt.scatter(neighborhood_matrix_knn[:, 0].T[mask], neighborhood_matrix_knn[:, 1].T[mask], color=color, label=str(color_cell_type_dict[i]), s=20, alpha=0.7)

            # plt.legend(fontsize=args.fontsize)
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
        
            # compute overlap of cell_type find found clusters
            computed_clusters = set(labels_clustering)

            with open('matches.txt', 'w') as matches_file:
                for i in computed_clusters:
                    mask_computed_clusters = labels_clustering == i
                    for j in range(len(cell_type_color_dict)):
                        mask_cell_type = labels_clustering_cell_type == j
                        
                        mask = mask_computed_clusters & mask_cell_type

                        number_of_matches = np.sum(mask)
                        matches_file.write('Computed cluster {} (size: {}) matching with cell type {} (size: {}) {} times. Rate (matches/computed_clusters): {}% Rate (matches/given_celltypes): {}%\n'.format(
                                        i, np.sum(mask_computed_clusters), color_cell_type_dict[j], np.sum(mask_cell_type), number_of_matches, number_of_matches/np.sum(mask_computed_clusters), number_of_matches/np.sum(mask_cell_type)))

                    matches_file.write('\n')
            
        if args.cell_coloring_batch:
            if len(colors) < len(cell_type_color_dict_batch):
                log.error('The chosen colormap offers too less values for the number of clusters.')
                exit(1)
            labels_clustering_cell_type_batch = []
            for cell_name in matrices_list:
                # if cell_name in cell_name_cell_type_dict:
                #     if cell_name_cell_type_dict[cell_name] in cell_type_color_dict:
                labels_clustering_cell_type_batch.append(cell_type_color_dict_batch[cell_name_cell_type_dict_batch[cell_name]])

            labels_clustering_cell_type_batch = np.array(labels_clustering_cell_type_batch)

            log.debug('labels_clustering_cell_type: {}'.format(len(labels_clustering_cell_type_batch)))
            log.debug('matrices_list: {}'.format(len(matrices_list)))


            plt.figure(figsize=(args.figuresize[0], args.figuresize[1]))
            for i, color in enumerate(colors[:len(cell_type_color_dict_batch)]):
                mask = labels_clustering_cell_type_batch == i
                log.debug('plot cluster: {} {}'.format(color_cell_type_dict_batch[i], np.sum(mask)))
                plt.scatter(neighborhood_matrix_knn[:, 0].T[mask], neighborhood_matrix_knn[:, 1].T[mask], color=color, label=str(color_cell_type_dict_batch[i]), s=20, alpha=0.7)

            # plt.legend(fontsize=args.fontsize)
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
        # plt.figure(figsize=(args.figuresize[0], args.figuresize[1]))
        # for i, color in enumerate(colors[:args.numberOfClusters]):
        #     mask = labels_clustering == i
        #     plt.scatter(neighborhood_matrix_knn[:, 1].T[mask], neighborhood_matrix_knn[:, 2].T[mask], color=color, label=str(i), s=50, alpha=0.7)
       
       
        ###### PCA 2 - 3 plots START
        # if args.cell_coloring:

        #     for i, color in enumerate(colors[:len(cell_type_color_dict)]):
        #         mask = labels_clustering_cell_type == i
        #         plt.scatter(neighborhood_matrix_knn[:, 1].T[mask], neighborhood_matrix_knn[:, 2].T[mask], color=color, label=str(color_cell_type_dict[i]), s=20, alpha=0.7)
        #     plt.legend(fontsize=args.fontsize)
        #     plt.xticks([])
        #     plt.yticks([])
        #     plt.xlabel('PC2', fontsize=args.fontsize)
        #     plt.ylabel('PC3', fontsize=args.fontsize)
        #     scatter_plot_name = '.'.join(args.createScatterPlot.split('.')[:-1]) + '_pc2_pc3_cell_color.' + args.createScatterPlot.split('.')[-1]
        #     plt.savefig(scatter_plot_name, dpi=args.dpi)
        #     plt.close()
        
        # plt.figure(figsize=(args.figuresize[0], args.figuresize[1]))
        # for i, color in enumerate(colors[:args.numberOfClusters]):
        #     mask = labels_clustering == i
        #     plt.scatter(neighborhood_matrix_knn[:, 1].T[mask], neighborhood_matrix_knn[:, 2].T[mask], color=color, label=str(i), s=20, alpha=0.7)
        # plt.legend(fontsize=args.fontsize)
        # plt.xticks([])
        # plt.yticks([])
        # plt.xlabel('PC2', fontsize=args.fontsize)
        # plt.ylabel('PC3', fontsize=args.fontsize)
        # scatter_plot_name = '.'.join(args.createScatterPlot.split('.')[:-1]) + '_pc2_pc3.' + args.createScatterPlot.split('.')[-1]
        # plt.savefig(scatter_plot_name, dpi=args.dpi)
        # plt.close()
        ###### PCA 2 - 3 plots END



        # plt.figure(figsize=(args.figuresize[0], args.figuresize[1]))
        # # for i, color in enumerate(colors[:args.numberOfClusters]):
        # #     mask = labels_clustering == i
        # #     plt.scatter(neighborhood_matrix_knn[:, 1].T[mask], neighborhood_matrix_knn[:, 2].T[mask], color=color, label=str(i), s=50, alpha=0.7)
        # if args.cell_coloring:

        #     for i, color in enumerate(colors[:len(cell_type_color_dict)]):
        #         mask = labels_clustering_cell_type == i
        #         plt.scatter(neighborhood_matrix_knn[:, 2].T[mask], neighborhood_matrix_knn[:, 3].T[mask], color=color, label=str(color_cell_type_dict[i]), s=20, alpha=0.7)
        # else:
        #     for i, color in enumerate(colors[:args.numberOfClusters]):
        #         mask = labels_clustering == i
        #         plt.scatter(neighborhood_matrix_knn[:, 2].T[mask], neighborhood_matrix_knn[:, 3].T[mask], color=color, label=str(i), s=20, alpha=0.7)
        # plt.legend(fontsize=args.fontsize)
        # plt.xticks([])
        # plt.yticks([])
        # plt.xlabel('PC2', fontsize=args.fontsize)
        # plt.ylabel('PC3', fontsize=args.fontsize)
        # scatter_plot_name = '.'.join(args.createScatterPlot.split('.')[:-1]) + '_pc3_pc4.' + args.createScatterPlot.split('.')[-1]
        # plt.savefig(scatter_plot_name, dpi=args.dpi)
        # plt.close()

    matrices_cluster = list(zip(matrices_list, labels_clustering))
    np.savetxt(args.outFileName, matrices_cluster, fmt="%s")
