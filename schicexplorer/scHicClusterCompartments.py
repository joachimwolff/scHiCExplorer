
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

from hicmatrix import HiCMatrix as hm
from schicexplorer.utilities import opener
from hicmatrix.lib import MatrixFileHandler
from hicexplorer.utilities import obs_exp_matrix_lieberman, obs_exp_matrix_norm, convertInfsToZeros_ArrayFloat
from hicexplorer.hicPCA import correlateEigenvectorWithGeneTrack, correlateEigenvectorWithHistonMarkTrack
from hicexplorer.utilities import convertNansToZeros, convertInfsToZeros
from sklearn.cluster import KMeans, SpectralClustering

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.cluster import SpectralClustering

from scipy.sparse import csr_matrix, lil_matrix
from scipy import linalg
from scipy.stats import pearsonr
import numpy as np
import pyBigWig

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
    parserRequired.add_argument('--distance', '-d',
                                help='Distance which should be considered as short range. Default 2MB.',
                                default=2000000,
                                type=int)
    parserRequired.add_argument('--chromosomes',
                                help='List of chromosomes to be included in the '
                                'correlation.',
                                default=None,
                                nargs='+')
    parserRequired.add_argument('--norm',
                                help='Different obs-exp normalization as used by '
                                'Homer software.',
                                action='store_true')
    parserRequired.add_argument('--binarization',
                                help='Set all positive values of eigenvetor to 1 and all negative ones to 0.',
                                action='store_true')
    parserRequired.add_argument('--extraTrack',
                                help='Either a gene track or a histon mark coverage'
                                ' file(preferably a broad mark) is needed to decide'
                                ' if the values of the eigenvector need a sign flip'
                                ' or not.',
                                default=None)
    parserRequired.add_argument('--histonMarkType',
                                help='set it to active or inactive. This is only '
                                'necessary if a histon mark coverage file is given '
                                'as an extraTrack.',
                                default='active')
    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the resulting clusters',
                                required=True,
                                default='clusters.txt')
    parserRequired.add_argument('--clusterMethod', '-cm',
                                help='Algorithm to cluster the Hi-C matrices',
                                choices=['spectral', 'kmeans'],
                                default='spectral')
    parserRequired.add_argument('--threads', '-t',
                                help='Number of threads. Using the python multiprocessing module.',
                                required=False,
                                default=4,
                                type=int)

    return parser


def open_and_store_matrix(pMatrixName, pMatricesList, pIndex, pXDimension, pChromosomes, pNorm, pExtraTrack, pHistonMarkType, pBinarization, pQueue):
    compartments_matrix = None
    for i, matrix in enumerate(pMatricesList):
        ma = hm.hiCMatrix(pMatrixName + '::' + matrix)
        ma.maskBins(ma.nan_bins)
        k = 1
        if pChromosomes:
            ma.keepOnlyTheseChr(pChromosomes)

        vecs_list = []
        chrom_list = []
        start_list = []
        end_list = []
        # PCA is computed per chromosome
        length_chromosome = 0
        chromosome_count = len(ma.getChrNames())

        for chrname in ma.getChrNames():
            chr_range = ma.getChrBinRange(chrname)
            length_chromosome += chr_range[1] - chr_range[0]
        if pExtraTrack and (pExtraTrack.endswith('.bw') or pExtraTrack.endswith('.bigwig')):
            bwTrack = pyBigWig.open(pExtraTrack, 'r')
        for chrname in ma.getChrNames():
            chr_range = ma.getChrBinRange(chrname)

            submatrix = ma.matrix[chr_range[0]:chr_range[1],
                                  chr_range[0]:chr_range[1]]
            if pNorm:
                obs_exp_matrix_ = obs_exp_matrix_norm(submatrix)

            else:
                obs_exp_matrix_ = obs_exp_matrix_lieberman(submatrix,
                                                           length_chromosome,
                                                           chromosome_count)
            obs_exp_matrix_ = convertNansToZeros(csr_matrix(obs_exp_matrix_)).todense()
            obs_exp_matrix_ = convertInfsToZeros(csr_matrix(obs_exp_matrix_)).todense()

            pearson_correlation_matrix = np.corrcoef(obs_exp_matrix_)
            pearson_correlation_matrix = convertNansToZeros(csr_matrix(pearson_correlation_matrix)).todense()
            pearson_correlation_matrix = convertInfsToZeros(csr_matrix(pearson_correlation_matrix)).todense()

            corrmatrix = np.cov(pearson_correlation_matrix)
            corrmatrix = convertNansToZeros(csr_matrix(corrmatrix)).todense()
            corrmatrix = convertInfsToZeros(csr_matrix(corrmatrix)).todense()
            evals, eigs = linalg.eig(corrmatrix)

            chrom, start, end, _ = zip(*ma.cut_intervals[chr_range[0]:chr_range[1]])

            chrom_list += chrom
            start_list += start
            end_list += end
            if pExtraTrack and (pExtraTrack.endswith('.bw') or pExtraTrack.endswith('.bigwig')):
                assert(len(end) == len(start))
                correlateEigenvectorWithHistonMarkTrack(eigs[:, :k].transpose(),
                                                        bwTrack, chrname, start,
                                                        end, pExtraTrack,
                                                        pHistonMarkType)

            vecs_list += eigs[:, :k].tolist()

        if compartments_matrix is None:
            compartments_matrix = np.zeros([pXDimension, len(np.array(vecs_list).flatten())], dtype=np.float)

        eigenvector = np.real(np.array(vecs_list).flatten())
        mask = np.isnan(eigenvector)
        if len(mask) > 0:
            eigenvector[mask] = 0
        mask = np.isinf(eigenvector)
        if len(mask) > 0:
            eigenvector[mask] = 0

        if pBinarization:
            mask = eigenvector <= 0
            eigenvector[mask] = -1
            mask = eigenvector > 0
            eigenvector[mask] = 1

        compartments_matrix[pIndex + i, :] = eigenvector

    pQueue.put(compartments_matrix)

    return


def main(args=None):

    args = parse_arguments().parse_args(args)

    create_or_load_matrix = False

    matrices_name = args.matrix
    threads = args.threads
    matrices_list = cooler.fileops.list_coolers(matrices_name)
    if threads > len(matrices_list):
        threads = len(matrices_list)
    compartments_matrix = None

    all_data_collected = False
    thread_done = [False] * threads
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
            pMatrixName=matrices_name,
            pMatricesList=matrices_name_list,
            pIndex=length_index[i],
            pXDimension=len(matrices_list),
            pChromosomes=args.chromosomes,
            pNorm=args.norm,
            pExtraTrack=args.extraTrack,
            pHistonMarkType=args.histonMarkType,
            pBinarization=args.binarization,
            pQueue=queue[i]
        )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                compartments_worker = queue[i].get()
                if compartments_matrix is None:
                    compartments_matrix = compartments_worker
                else:
                    compartments_matrix += compartments_worker

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
        spectral_clustering = SpectralClustering(n_clusters=args.numberOfClusters, n_jobs=args.threads)
        labels_clustering = spectral_clustering.fit_predict(compartments_matrix)
    elif args.clusterMethod == 'kmeans':
        kmeans_object = KMeans(n_clusters=args.numberOfClusters, random_state=0, n_jobs=args.threads, precompute_distances=True)
        labels_clustering = kmeans_object.fit_predict(compartments_matrix)

    matrices_cluster = list(zip(matrices_list, labels_clustering))
    np.savetxt(args.outFileName, matrices_cluster, fmt="%s")
