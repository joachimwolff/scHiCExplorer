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
    parserRequired.add_argument('--outFileName', '-o',
                                help='File name of the consensus mcool matrix.',
                                required=True)

    parserRequired.add_argument('--threads',
                           help='Number of threads. Using the python multiprocessing module.',
                           required=False,
                           default=4,
                           type=int)

    return parser



def compute_consensus_matrix(pMatrixName, pMatricesList, pAppend, pCounter, pQueue):
    read_coverage = []
    sparsity = []
    read_distribution = []
    max_length = 0
    consensus_matrix = None
    for i, matrix in enumerate(pMatricesList):
        matrixFileHandlerInput = MatrixFileHandler(pFileType='cool', pMatrixFile=pMatrixName+ '::' +matrix)
        _matrix, cut_intervals, nan_bins, \
                distance_counts, correction_factors = matrixFileHandlerInput.load()
    
        if consensus_matrix is None:
            consensus_matrix = _matrix
        else:
            consensus_matrix += _matrix
    
    hic2CoolVersion = matrixFileHandlerInput.matrixFile.hic2cool_version
    matrixFileHandlerOutput = MatrixFileHandler(pFileType='cool', pAppend=pAppend, pEnforceInteger=False, pFileWasH5=False, pHic2CoolVersion=hic2CoolVersion)

    matrixFileHandlerOutput.set_matrix_variables(consensus_matrix, cut_intervals, nan_bins,
                                                    correction_factors, distance_counts)
    # matrixFileHandlerOutput.save(
    #     pOutFileName + '::/' + pConsensusMatrixName, pSymmetric=True, pApplyCorrection=False)
    pQueue.put([matrixFileHandlerOutput, pCounter])

def main(args=None):

    args = parse_arguments().parse_args(args)

    matrices_name = args.matrix
    matrices_list = cooler.fileops.list_coolers(matrices_name)
    clusters = {}
    with open(args.clusters, 'r') as cluster_file:

        for i, line in enumerate(cluster_file.readlines()):
            line = line.strip()
            file_path, cluster = line.split(' ')
            
            if int(cluster) in clusters:
                clusters[int(cluster)].append(file_path)
            else:
                clusters[int(cluster)] = [file_path]

    cluster_list = []
    for key in clusters:
        cluster_list.append(clusters[key])

    threads = args.threads

    process = [None] * args.threads
    all_data_processed = False
    # hic_matrix = coo_matrix((matrix_size, matrix_size), dtype='uint32')
    queue = [None] * args.threads

    all_threads_done = False
    thread_done = [False] * args.threads
    count_output = 0
    count_call_of_read_input = 0
    computed_pairs = 0

    matrixFileHandlerObjects_list = [None] * len(cluster_list)

    while not all_data_processed or not all_threads_done:

        for i in range(args.threads):
            if queue[i] is None and not all_data_processed:
                if count_call_of_read_input < len(cluster_list):
                    queue[i] = Queue()
                    process[i] = Process(target=compute_consensus_matrix, kwargs=dict(
                                        pMatrixName = matrices_name,
                                        pMatricesList= cluster_list[count_call_of_read_input], 
                                        # pOutFileName=args.outFileName,
                                        # pConsensusMatrixName = 'consensus_matrix_cluster' + str(count_call_of_read_input), 
                                        pAppend = i > 0, 
                                        pCounter = count_call_of_read_input,
                                        pQueue=queue[i]
                        )
                    )

                    process[i].start()
                    thread_done[i] = False
                else:
                    all_data_processed = True
                    thread_done[i] = True
                count_call_of_read_input += 1
            elif queue[i] is not None and not queue[i].empty():
                matrixFileHandlerObjects, counter_ = queue[i].get()
                matrixFileHandlerObjects_list[counter_] = matrixFileHandlerObjects
                queue[i] = None
                process[i].join()
                process[i].terminate()
                process[i] = None
                thread_done[i] = True
            elif all_data_processed and queue[i] is None:
                thread_done[i] = True
            else:
                time.sleep(1)

        if all_data_processed:
            all_threads_done = True
            for thread in thread_done:
                if not thread:
                    all_threads_done = False


    for i, matrixFileHandler in enumerate(matrixFileHandlerObjects_list):
        matrixFileHandler.save(
            args.outFileName + '::/' + 'consensus_matrix_cluster_' + str(i), pSymmetric=True, pApplyCorrection=False)