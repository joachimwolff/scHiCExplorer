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

from hicmatrix import HiCMatrix
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
from mpl_toolkits.axes_grid1 import AxesGrid

def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        add_help=False
    )

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--matrix', '-m',
                                help='The single cell Hi-C interaction matrices to investigate for QC. Needs to be in mcool format',
                                metavar='mcool scHi-C matrix',
                                required=True)
    parserRequired.add_argument('--clusters', '-c',
                                help='Text file which contains per matrix the assoziated cluster.',
                                metavar='cluster file',
                                required=True)
    parserRequired.add_argument('--maximalDistance',
                           help='maximal distance in bases',
                           required=False,
                           default=50000000,
                           type=int)

    parserRequired.add_argument('--threads',
                           help='Number of threads. Using the python multiprocessing module.',
                           required=False,
                           default=4,
                           type=int)

    return parser



def compute_read_distribution(pMatrixName, pMatricesList, pIndex, pXDimension, pQueue):
    read_coverage = []
    sparsity = []
    read_distribution = []
    max_length = 0
    for i, matrix in enumerate(pMatricesList):
        matrixFileHandlerInput = MatrixFileHandler(pFileType='cool', pMatrixFile=pMatrixName+ '::' +matrix)
        _matrix, _, _, _, _ = matrixFileHandlerInput.load()
        read_distribution_ = {}
        # if neighborhood_matrix is None:
        #     neighborhood_matrix = csr_matrix((pXDimension, _matrix.shape[0] * _matrix.shape[1]), dtype=np.float)


        instances, features = _matrix.nonzero()

        relative_distance = np.absolute(instances - features)
        if max_length < max(relative_distance):
            max_length = max(relative_distance)
        
        for i, relative_distance_ in enumerate(relative_distance):
            if relative_distance_ in read_distribution_:
                read_distribution_[relative_distance_] += _matrix.data[i]
            else:
                read_distribution_[relative_distance_] = _matrix.data[i]

        read_distribution.append(read_distribution_)
    
    pQueue.put([read_distribution, max_length])

def main(args=None):

    args = parse_arguments().parse_args(args)

    create_or_load_matrix = False

    # if args.createMatrix:
        
    matrices_name = args.matrix
    threads = args.threads
    matrices_list = cooler.fileops.list_coolers(matrices_name)
    read_coverage = [None] * threads
    sparsity = [None] * threads


    all_data_collected = False
    thread_done = [False] * threads
    log.debug('matrix read, starting processing')
    length_index = [None] * threads
    length_index[0] = 0
    matricesPerThread = len(matrices_list) // threads
    queue = [None] * threads
    process = [None] * threads
    max_length = 0
    for i in range(threads):

        if i < threads - 1:
            matrices_name_list = matrices_list[i * matricesPerThread:(i + 1) * matricesPerThread]
            length_index[i + 1] = length_index[i] + len(matrices_name_list)
        else:
            matrices_name_list = matrices_list[i * matricesPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=compute_read_distribution, kwargs=dict(
                            pMatrixName = matrices_name,
                            pMatricesList= matrices_name_list, 
                            pIndex = length_index[i], 
                            pXDimension=len(matrices_list),
                            pQueue=queue[i]
            )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                read_coverage[i], max_length_ = queue[i].get()
                if max_length < max_length_:
                    max_length = max_length_
                # read_coverage[i] = worker_result
                # sparsity[i] = worker_result[1]

                # if neighborhood_matrix is None:
                #     neighborhood_matrix = csr_matrix_worker
                # else:
                #     neighborhood_matrix += csr_matrix_worker

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

    read_coverage = np.array([item for sublist in read_coverage for item in sublist])

    clusters = {}
    with open(args.clusters, 'r') as cluster_file:

        for i, line in enumerate(cluster_file.readlines()):
            line = line.strip()
            line_ = line.split(' ')[1]
            
            if int(line_) in clusters:
                clusters[int(line_)].append(read_coverage[i])
            else:
                clusters[int(line_)] = [read_coverage[i]]

    cluster_to_plot = []    
    # max_length = 0
    # for key in clusters:
    #     if max_length < max(clusters[key]):
    #         max_length = max(clusters[key])

    w=10
    h=100
    # fig=plt.figure(figsize=(20, 10))
    columns = 12
    rows = 1
    # for i, key in enumerate(background_model_data):
    # #     img = np.random.randint(10, size=(h,w))
    #     fig.add_subplot(rows, columns, i+1)
    # #     plt.imshow(img)
    #     plt.hist(x=background_model_data[key], bins=10, log=True)
    # plt.savefig('histograms_nolog.png', dpi=300)
    # fig = plt.figure()

    # grid = AxesGrid(fig, 111,
    #                 nrows_ncols=(1, 12),
    #                 axes_pad=0.005,
    #                 share_all=True,
    #                 label_mode="L",
    #                 cbar_location="right",
    #                 cbar_mode="single",
    #                 )

    # for val, ax in zip(vals,grid):

   
    binSize = 1000000

    if (args.maximalDistance // binSize) < max_length:

        max_length = args.maximalDistance // binSize
    # fig, axes = plt.subplots(rows, columns, sharey=True)
    
    clusters_list = []
    cluster_size = []
    for i, key in enumerate(clusters):
        cluster_to_plot = []
        
        # log.debug('cluster length {}'.format(len(clusters[key])))
        for cluster_item in clusters[key]:
            cluster_array = np.zeros(max_length)
            for distance in cluster_item.keys():
                if distance < max_length:
                    cluster_array[distance] = cluster_item[distance]
            
            # cluster_array.append(np.array(list(cluster_item.values()), dtype=np.int))
            cluster_to_plot.append(cluster_array)
        # im = grid[i].imshow(np.log10(np.array(cluster_to_plot).T), cmap='hot', interpolation='nearest')
        clusters_list.append(cluster_to_plot)
        cluster_size.append(len(cluster_to_plot))

    cluster_size = np.array(cluster_size)
    log.debug('cluster size {}'.format(cluster_size))
    cluster_size = (cluster_size-min(cluster_size))/(max(cluster_size)-min(cluster_size))
    log.debug('cluster size {}'.format(cluster_size))
    cluster_size = list(cluster_size)
    cluster_size.append(0.08)
    f,axes = plt.subplots(rows, columns, 
            gridspec_kw={'width_ratios':cluster_size})
    for i in range(1, len(axes)-1):
        axes[0].get_shared_y_axes().join(axes[0],axes[i])
    # ax1.get_shared_y_axes().join(ax2,ax3)
    for i, cluster in enumerate(clusters_list):
        axes[i].imshow(np.log10(np.array(cluster).T), cmap='hot', interpolation='nearest')
        # plt.invert_yaxis()
        # plt.xaxis.set_visible(False)
        # if i % 2:
        #     g.yaxis.set_visible(False)
        # # plt.axes.Axes.set_xscale(1, 'linear')
        # # plt.legend('Cluster {}'.format(i))

        # plt.yscale('log')
    # grid[0].yaxis.set_visible(True)

    # grid[-1].yaxis.set_visible(True)
    # grid.cbar_axes[0].colorbar(im)

    # for cax in grid.cbar_axes:
    #     cax.toggle_label(False)
    # log.debug('cluster_array {}'.format(cluster_array))
    # for i in cluster_to_plot:
    # log.debug('clusters {}'.format(cluster_to_plot))

    # plt.imshow(np.log10(cluster_to_plot).T, cmap='hot', interpolation='nearest')
    plt.savefig('density_plot.png', dpi=300)