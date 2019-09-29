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
    parserRequired.add_argument('--loadFromNpz', '-npz',
                                help='Load distance data from npz file.',
                                # metavar='Create npz matrix or load it.',
                                # required=True,
                                action='store_false')
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

def getIndex(pValue, pDistances):

    # distances 10, 100, 1000, 10000
    # for 10: 0,1, 2,3 ,4,5,6,7,8,9 --> same index pos
    # for 100: merge 10 - 19 to one, 20-29 etc --> index pis 10, 11, 12
    # for 1000: merge 100-199 to one etc --> index pos 20, 21
    if pValue < pDistances[0]:
        return pValue
    for i, distance in enumerate(pDistances):
        if pValue <= distance:
            rest_factor = pValue // (pDistances[i-1])
            ten_factor = np.floor(np.log10(pValue)) * 10
            # log.debug('pValue {} {}'.format(np.log10(pValue)))

            # log.debug('pValue {} int(ten_factor + rest_factor) {} {}'.format(pValue, ten_factor, rest_factor))
            return int(ten_factor + rest_factor)
    
    return -1

def compute_read_distribution(pMatrixName, pMatricesList, pIndex, pXDimension, pSteps, pQueue):
    read_coverage = []
    sparsity = []
    read_distribution = []
    max_length = 0
    for i, matrix in enumerate(pMatricesList):
        matrixFileHandlerInput = MatrixFileHandler(pFileType='cool', pMatrixFile=pMatrixName+ '::' +matrix)
        _matrix, _, _, _, _ = matrixFileHandlerInput.load()
        
        # if neighborhood_matrix is None:
        #     neighborhood_matrix = csr_matrix((pXDimension, _matrix.shape[0] * _matrix.shape[1]), dtype=np.float)


        instances, features = _matrix.nonzero()

        relative_distance = np.absolute(instances - features)
        if max_length < max(relative_distance):
            max_length = max(relative_distance)
        read_distribution_ = np.zeros(max_length+1)
        for j, relative_distance_ in enumerate(relative_distance):
            read_distribution_[relative_distance_] += _matrix.data[j]
            # if relative_distance_ == 0:
            #     index = 0
            # else:
            #     index = getIndex(relative_distance_, pSteps)
            # if index == -1:
            #     read_distribution_[-1] += _matrix.data[j]
            #     continue
            #     # cluster_array[index] = cluster_item[distance]
            # if max_index < index:
            #     max_index = index
            # if index in read_distribution_:
            #     read_distribution_[index] += _matrix.data[j]
            # else:
            #     read_distribution_[index] = _matrix.data[j]

        read_distribution.append(read_distribution_)
    
    pQueue.put([read_distribution, max_length])

def main(args=None):

    args = parse_arguments().parse_args(args)

    create_or_load_matrix = False
    binSize = 10000
    max_length = 100000000 # 100 mb
    max_index = 0
    steps = np.array([100000, 1000000, 10000000, 100000000 ])
    steps = steps / binSize
    if args.loadFromNpz:
        
        matrices_name = args.matrix
        threads = args.threads
        matrices_list = cooler.fileops.list_coolers(matrices_name)
        log.debug('matrices_list {}'.format(matrices_list[:10]))
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
                                pSteps=steps,
                                pQueue=queue[i]
                )
            )

            process[i].start()

        while not all_data_collected:
            for i in range(threads):
                if queue[i] is not None and not queue[i].empty():
                    read_coverage[i], max_index_ = queue[i].get()
                    if max_index < max_index_:
                        max_index = max_index_
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

        log.debug('matrices processed')
        read_distributions = []
        for thread_data in read_coverage:
            for matrix_data in thread_data:
                # matrix_data_array = [0] * (max_index + 1)
                # for key, value in matrix_data.items():
                #     matrix_data_array[key] = value
                read_distributions.append(matrix_data)
        
        read_distributions = np.array(read_distributions)
        # read_coverage = np.array([item for sublist in read_coverage for item in sublist])
        log.debug('one array computation done')
        np.save('readcoverage_data.npz', np.array(read_distributions))

    else:
        log.debug('read file')
        read_distributions = np.load('readcoverage_data.npz.npy', allow_pickle=True)

    log.debug('Assoziating clusters')
    clusters = {}
    with open(args.clusters, 'r') as cluster_file:

        for i, line in enumerate(cluster_file.readlines()):
            line = line.strip()
            line_ = line.split(' ')[1]
            
            if int(line_) in clusters:
                clusters[int(line_)].append(read_distributions[i])
            else:
                clusters[int(line_)] = [read_distributions[i]]

    cluster_to_plot = []    

    log.debug('len(read_distributions) {}'.format(len(read_distributions)))
    log.debug('len(read_distributions[0]) {}'.format(len(read_distributions[0])))

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


   
    length_plot = 10 * len(steps)
    if (args.maximalDistance // binSize) < max_length:

        max_length = args.maximalDistance // binSize
    # fig, axes = plt.subplots(rows, columns, sharey=True)
    log.debug('Pooling data')
    clusters_list = []
    cluster_size = []
    for i, key in enumerate(clusters):
        cluster_to_plot = []
        
        # log.debug('cluster length {}'.format(len(clusters[key])))
        for cluster_item in clusters[key]:
            # cluster_array = np.zeros(max_length)
            # log.debug('len(clusters[key]) {}'.format(len(clusters[key])))
            # for distance in cluster_item.keys():
               


            #     if distance < max_length:
            #     # log_value_distance = np.floor(np.log10(distance))
            #     # if log_value_distance < 2:
            #         cluster_array[distance] = cluster_item[distance]
              
            
            # cluster_array.append(np.array(list(cluster_item.values()), dtype=np.int))
            cluster_to_plot.append(cluster_item)
        # im = grid[i].imshow(np.log10(np.array(cluster_to_plot).T), cmap='hot', interpolation='nearest')
        clusters_list.append(cluster_to_plot)
        cluster_size.append(len(cluster_to_plot))

        # np.savez('cluster_data.npz', np.array(clusters_list))
    # else:
    #     clusters_list = np.load(args.matrix)
    cluster_size = np.array(cluster_size)
    # print(len(cluster_size))
    cluster_size = (1.0 - 0.1) * (cluster_size-min(cluster_size))/(max(cluster_size)-min(cluster_size)) + (0.1)
    # print(len(cluster_size))
    cluster_size = list(cluster_size)
    log.debug('Plot data')

    # print(len(cluster_size))
    f,axes = plt.subplots(rows, columns, 
            gridspec_kw={'width_ratios':cluster_size})
    for i, cluster in enumerate(clusters_list):
        axes[i].imshow(np.array(cluster).T, interpolation='nearest')
        log.debug('np.array(cluster[0]).T {}'.format(np.array(cluster[0])))
        axes[i].set_yscale('log')
        axes[i].invert_yaxis()
        axes[i].get_xaxis().set_ticks([])
        if i > 0:
            axes[i].yaxis.set_visible(False)
        axes[i].set_xlabel(str(i))
    plt.savefig('density_plot.png', dpi=300)

    plt.close()

    colors = ['silver', 'gray', 'black', 'red', 'maroon', 'yellow', 'olive', 'lime', 'green', 'aqua', 'teal', 'blue', 'navy', 'fuchsia', 'purple']

    # for i, cluster in enumerate(clusters_list):
    #     for matrix_data in cluster:
    #         plt.plot(list(range(len(matrix_data))), matrix_data, c=colors[i%len(colors)], alpha=0.5)
        
    # plt.yscale('log')

    # plt.savefig('distance_coverage.png', dpi=300)
    