
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
logging.getLogger('cooler').setLevel(logging.WARNING)
logging.getLogger('hicmatrix').setLevel(logging.WARNING)


log = logging.getLogger(__name__)

import cooler

from hicmatrix import HiCMatrix
from schicexplorer.utilities import opener
from hicmatrix.lib import MatrixFileHandler
from hicexplorer.utilities import obs_exp_matrix_lieberman, obs_exp_matrix_norm

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.cluster import SpectralClustering


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
    parserRequired.add_argument('--matrixNpz', '-npz',
                                help='The single cell Hi-C interaction matrices to cluster. Needs to be in mcool format',
                                metavar='mcool scHi-C matrix',
                                required=False)
    parserRequired.add_argument('--createMatrix', '-cm',
                                help='If set to, the matrix for the clustering is created out of the single cell mcool matrix. If not, the binary npz matrix of a former creation is loaded.',
                                # metavar='Create npz matrix or load it.',
                                # required=True,
                                action='store_true')
    parserRequired.add_argument('--numberOfClusters', '-c',
                           help='Number of to be computed clusters',
                           required=False,
                           default=12,
                           type=int)
    parserRequired.add_argument('--distance', '-d',
                           help='Distance which should be considered as short range. Default 2MB.',
                           default=2000000,
                           type=int)

    parserRequired.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module.',
                           required=False,
                           default=4,
                           type=int)

    return parser



def open_and_store_matrix(pMatrixName, pMatricesList, pIndex, pXDimension, pNorm, pQueue):
    compartments_matrix = None
    for i, matrix in enumerate(pMatricesList):
        # matrixFileHandlerInput = MatrixFileHandler(pFileType='cool', pMatrixFile=pMatrixName+ '::' +matrix)
        # _matrix, _, _, _, _ = matrixFileHandlerInput.load()
        # chromosomes_list = cooler.Cooler(pMatrixName+ '::' +matrix).chromnames
        ma = hm.hiCMatrix(pMatrixName+ '::' +matrix)
        ma.maskBins(ma.nan_bins)
        k = 1
        if args.chromosomes:
            ma.keepOnlyTheseChr(args.chromosomes)

        vecs_list = []
        chrom_list = []
        start_list = []
        end_list = []
        # PCA is computed per chromosome
        length_chromosome = 0
        chromosome_count = len(ma.getChrNames())
        if args.pearsonMatrix:
            trasf_matrix_pearson = lil_matrix(ma.matrix.shape)

        if args.obsexpMatrix:
            trasf_matrix_obsexp = lil_matrix(ma.matrix.shape)

        for chrname in ma.getChrNames():
            chr_range = ma.getChrBinRange(chrname)
            length_chromosome += chr_range[1] - chr_range[0]
        if args.extraTrack and (args.extraTrack.endswith('.bw') or args.extraTrack.endswith('.bigwig')):
            bwTrack = pyBigWig.open(args.extraTrack, 'r')
        for chrname in ma.getChrNames():
            chr_range = ma.getChrBinRange(chrname)

            submatrix = ma.matrix[chr_range[0]:chr_range[1],
                                chr_range[0]:chr_range[1]]
            if args.norm:
                obs_exp_matrix_ = obs_exp_matrix_norm(submatrix)

            else:
                obs_exp_matrix_ = obs_exp_matrix_lieberman(submatrix,
                                                        length_chromosome,
                                                        chromosome_count)
            obs_exp_matrix_ = convertNansToZeros(csr_matrix(obs_exp_matrix_)).todense()
            obs_exp_matrix_ = convertInfsToZeros(csr_matrix(obs_exp_matrix_)).todense()

            if args.obsexpMatrix:
                trasf_matrix_obsexp[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(obs_exp_matrix_)

            pearson_correlation_matrix = np.corrcoef(obs_exp_matrix_)
            pearson_correlation_matrix = convertNansToZeros(csr_matrix(pearson_correlation_matrix)).todense()
            pearson_correlation_matrix = convertInfsToZeros(csr_matrix(pearson_correlation_matrix)).todense()

            if args.pearsonMatrix:
                trasf_matrix_pearson[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(pearson_correlation_matrix)

            corrmatrix = np.cov(pearson_correlation_matrix)
            corrmatrix = convertNansToZeros(csr_matrix(corrmatrix)).todense()
            corrmatrix = convertInfsToZeros(csr_matrix(corrmatrix)).todense()
            evals, eigs = linalg.eig(corrmatrix)
            k = args.numberOfEigenvectors

            chrom, start, end, _ = zip(*ma.cut_intervals[chr_range[0]:chr_range[1]])

            chrom_list += chrom
            start_list += start
            end_list += end
            if args.extraTrack and (args.extraTrack.endswith('.bw') or args.extraTrack.endswith('.bigwig')):
                assert(len(end) == len(start))
                correlateEigenvectorWithHistonMarkTrack(eigs[:, :k].transpose(),
                                                        bwTrack, chrname, start,
                                                        end, args.extraTrack,
                                                        args.histonMarkType)

            vecs_list += eigs[:, :k].tolist()


        if compartments_matrix is None:
            compartments_matrix = csr_matrix((pXDimension, len(chromosomes_list)), dtype=np.float)

        compartments_matrix[pIndex+i, :] = np.array(vecs_list)
        # if i % 20 == 0:
        #     log.debug('pIndex + i {} {}'.format(pIndex, i))
    
    pQueue.put(compartments_matrix)

    return

def main(args=None):

    args = parse_arguments().parse_args(args)

    create_or_load_matrix = False

    if args.createMatrix:
        
        matrices_name = args.matrix
        threads = args.threads
        matrices_list = cooler.fileops.list_coolers(matrices_name)
        svl_matrix = None

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
                                pQueue=queue[i]
                )
            )

            process[i].start()

        while not all_data_collected:
            for i in range(threads):
                if queue[i] is not None and not queue[i].empty():
                    csr_matrix_worker = queue[i].get()
                    if svl_matrix is None:
                        svl_matrix = csr_matrix_worker
                        log.debug('returned first csr i {}'.format(i))
                    else:
                        svl_matrix += csr_matrix_worker
                        log.debug('adding csr i {}'.format(i))

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
        scipy.sparse.save_npz(args.matrix + '_binary.npz', svl_matrix)
    else:
        log.debug('read npz')
        svl_matrix = scipy.sparse.load_npz(args.matrixNpz)
        matrices_list = cooler.fileops.list_coolers(args.matrix)


    
    spectral_clustering = SpectralClustering(n_clusters=args.numberOfClusters, n_jobs=args.threads)
    labels_clustering = spectral_clustering.fit_predict(svl_matrix)

    matrices_cluster = list(zip(matrices_list, labels_clustering))
    np.savetxt('matrices_cluster_svl.txt', matrices_cluster, fmt="%s")