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
    parserRequired.add_argument('--createMatrix', '-cm',
                                help='If set to, the matrix for the clustering is created out of the single cell mcool matrix. If not, the binary npz matrix of a former creation is loaded.',
                                # metavar='Create npz matrix or load it.',
                                # required=True,
                                action='store_true')
    parserRequired.add_argument('--outputMcool', '-o',
                                help='Mcool matrix which contains only the filtered matrices',
                                
                                default='clusters')

    parserRequired.add_argument('--threads',
                           help='Number of threads. Using the python multiprocessing module.',
                           required=False,
                           default=4,
                           type=int)

    return parser



def compute_read_coverage_sparsity(pMatrixName, pMatricesList, pIndex, pXDimension, pQueue):
    read_coverage = []
    sparsity = []
    for i, matrix in enumerate(pMatricesList):
        matrixFileHandlerInput = MatrixFileHandler(pFileType='cool', pMatrixFile=pMatrixName+ '::' +matrix)
        _matrix, _, _, _, _ = matrixFileHandlerInput.load()

        # if neighborhood_matrix is None:
        #     neighborhood_matrix = csr_matrix((pXDimension, _matrix.shape[0] * _matrix.shape[1]), dtype=np.float)

        read_coverage.append(_matrix.data.sum())
        sparsity.append(_matrix.nnz / (_matrix.shape[0] * _matrix.shape[1]))

        # instances, features = _matrix.nonzero()

        # instances *= _matrix.shape[1]
        # instances += features
        # features = None
        # neighborhood_matrix[pIndex+i, instances] = _matrix.data
    
    pQueue.put([read_coverage, sparsity])

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
    for i in range(threads):

        if i < threads - 1:
            matrices_name_list = matrices_list[i * matricesPerThread:(i + 1) * matricesPerThread]
            length_index[i + 1] = length_index[i] + len(matrices_name_list)
        else:
            matrices_name_list = matrices_list[i * matricesPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=compute_read_coverage_sparsity, kwargs=dict(
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
                worker_result = queue[i].get()
                read_coverage[i] = worker_result[0]
                sparsity[i] = worker_result[1]

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
    sparsity = np.array([item for sublist in sparsity for item in sublist])
   

    magnitudes = np.log10(read_coverage)
    average_magnitude = np.int(np.average(magnitudes))

    magnitudes_sparsity = np.log10(sparsity)
    sparsity_average_magnitude = np.int(np.average(magnitudes_sparsity))
    
    y = [0] * len(read_coverage)
    # for i in range(len(args.matrices)):
    #     y[i] = [i] * len(sparsity)
    plt.axvline(x=average_magnitude, color='r')
    plt.plot(np.log10(read_coverage), y,  'o', alpha=0.5, markersize=0.3)
    plt.savefig('read_coverage_distribution.png', dpi=300)
    plt.close()
    plt.plot(sparsity, y,  'o', alpha=0.5, markersize=0.3)
    # plt.plot(np.log10(read_coverage), y,  'o', alpha=0.5, markersize=0.3)

    plt.savefig('sparsity_distribution.png', dpi=300)
    plt.close()
    plt.axvline(x=sparsity_average_magnitude, color='r')

    plt.plot(np.log10(sparsity), y,  'o', alpha=0.5, markersize=0.3)
    # plt.plot(np.log10(read_coverage), y,  'o', alpha=0.5, markersize=0.3)

    plt.savefig('sparsity_log10_distribution.png', dpi=300)
    plt.close()

    # scipy.sparse.save_npz(args.matrix + 'binary.npz', neighborhood_matrix)
    # else:
    #     neighborhood_matrix = scipy.sparse.load_npz(args.matrix)

    mask_average_magnitude = magnitudes >= average_magnitude
    mask_average_sparsity = magnitudes_sparsity >= sparsity_average_magnitude
    mask = np.logical_or(mask_average_magnitude, mask_average_sparsity)
    # read_coverage_filtered = read_coverage[mask]
    matrices_list_filtered = np.array(matrices_list)[mask]
    log.debug('len(matrices_list_filtered) {}'.format(len(matrices_list_filtered)))
    log.debug('matrices_list_filtered {}'.format(matrices_list_filtered))
    np.savetxt('accepted_matrices.txt', matrices_list_filtered,  fmt="%s")
    np.savetxt('rejected_matrices.txt', np.array(matrices_list)[~mask], fmt="%s")

    for matrix in matrices_list_filtered:

        cooler.fileops.cp(args.matrix +'::' + matrix, args.outputMcool + '::' + matrix)