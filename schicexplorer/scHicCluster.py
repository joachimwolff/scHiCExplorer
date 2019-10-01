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

from sklearn.cluster import KMeans, SpectralClustering

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
    parserRequired.add_argument('--numberOfClusters', '-c',
                           help='Number of to be computed clusters',
                           required=False,
                           default=12,
                           type=int)
    parserRequired.add_argument('--numberOfNeighbors', '-n',
                           help='Number of neighbors of clustering',
                           required=False,
                           default=10,
                           type=int)

    parserRequired.add_argument('--chromosomes',
                           help='List of to be plotted chromosomes',
                           nargs='+')
    parserRequired.add_argument('--clusterMethod', '-cm',
                           help='Algorithm to cluster the Hi-C matrices',
                           choices=['spectral', 'kmeans'],
                           default='spectral')
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
        # matrixFileHandlerInput = MatrixFileHandler(pFileType='cool', pMatrixFile=pMatrixName+ '::' +matrix)
        # _matrix, _, _, _, _ = matrixFileHandlerInput.load()
        _matrix = hic_ma.matrix

        if neighborhood_matrix is None:
            neighborhood_matrix = csr_matrix((pXDimension, _matrix.shape[0] * _matrix.shape[1]), dtype=np.float)

        instances, features = _matrix.nonzero()

        instances *= _matrix.shape[1]
        instances += features
        features = None
        neighborhood_matrix[pIndex+i, instances] = _matrix.data
        # if i % 20 == 0:
        #     log.debug('pIndex + i {} {}'.format(pIndex, i))
    
    pQueue.put(neighborhood_matrix)

def main(args=None):

    args = parse_arguments().parse_args(args)

    create_or_load_matrix = False

    # if args.createMatrix:
        
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
                    # log.debug('returned first csr i {}'.format(i))
                else:
                    neighborhood_matrix += csr_matrix_worker
                    # log.debug('adding csr i {}'.format(i))

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

    if args.clusterMethod == 'spectral':
        spectralClustering_object = SpectralClustering(n_clusters=args.numberOfClusters, n_jobs=args.threads,
                                                        n_neighbors=args.numberOfNeighbors)

        labels_clustering = spectralClustering_object.fit_predict(neighborhood_matrix)
    elif args.clusterMethod == 'kmeans':
        log.debug('start kmeans')

        kmeans_object = KMeans(n_clusters=args.numberOfClusters, random_state=0, n_jobs=args.threads)

        labels_clustering = kmeans_object.fit_predict(neighborhood_matrix)

    matrices_cluster = list(zip(matrices_list, labels_clustering))
    np.savetxt(args.outFileName, matrices_cluster, fmt="%s")
