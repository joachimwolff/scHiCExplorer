
import argparse
import os
import gzip
import shutil
from multiprocessing import Process, Queue
import time

import logging
log = logging.getLogger(__name__)

import cooler

from hicmatrix import HiCMatrix as hm

import numpy as np
from scipy.sparse import csr_matrix
from sklearn.cluster import SpectralClustering, KMeans


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
    parserRequired.add_argument('--distanceShortRange', '-ds',
                                help='Distance which should be considered as short range. Default 2MB.',
                                default=2000000,
                                type=int)
    parserRequired.add_argument('--distanceLongRange', '-dl',
                                help='Distance which should be considered as short range. Default 12MB.',
                                default=12000000,
                                type=int)
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


def create_svl_data(pMatrixName, pMatricesList, pIndex, pXDimension, pDistanceMin, pDistanceMax, pQueue):
    svl_matrix = None
    for i, matrix in enumerate(pMatricesList):
        chromosomes_list = cooler.Cooler(pMatrixName + '::' + matrix).chromnames
        svl_relations = []

        for chromosome in chromosomes_list:
            hic_matrix_obj = hm.hiCMatrix(
                pMatrixFile=pMatrixName + '::' + matrix, pChrnameList=[chromosome])
            if hic_matrix_obj.matrix.shape[0] < 5:
                svl_relations.append(0)
                continue
            min_distance = pDistanceMin / hic_matrix_obj.getBinSize()
            max_distance = pDistanceMax / hic_matrix_obj.getBinSize()

            hic_matrix = hic_matrix_obj.matrix

            instances, features = hic_matrix.nonzero()
            distances = np.absolute(instances - features)
            mask = distances <= min_distance
            mask_mitotic_0 = distances > min_distance
            mask_mitotic_1 = max_distance <= distances

            mask_mitotic = np.logical_and(mask_mitotic_0, mask_mitotic_1)

            sum_smaller_max_distance = np.sum(hic_matrix.data[mask])
            sum_greater_max_distance = np.sum(hic_matrix.data[mask_mitotic])
            svl_relation = sum_smaller_max_distance / sum_greater_max_distance
            if np.isinf(svl_relation) or np.isnan(svl_relation):
                svl_relations.append(0)
                continue
            svl_relations.append(svl_relation)

        if svl_matrix is None:
            svl_matrix = csr_matrix((pXDimension, len(chromosomes_list)), dtype=np.float)

        svl_matrix[pIndex + i, :] = np.array(svl_relations)

    pQueue.put(svl_matrix)

    return


def main(args=None):

    args = parse_arguments().parse_args(args)

    matrices_name = args.matrix
    threads = args.threads
    matrices_list = cooler.fileops.list_coolers(matrices_name)
    svl_matrix = None

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
        process[i] = Process(target=create_svl_data, kwargs=dict(
            pMatrixName=matrices_name,
            pMatricesList=matrices_name_list,
            pIndex=length_index[i],
            pXDimension=len(matrices_list),
            pDistanceMin=args.distanceShortRange,
            pDistanceMax=args.distanceLongRange,
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
                else:
                    svl_matrix += csr_matrix_worker

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
        log.debug('SVL matrix created, start clustering spectral')
        spectral_clustering = SpectralClustering(n_clusters=args.numberOfClusters, affinity='nearest_neighbors', n_jobs=args.threads)
        labels_clustering = spectral_clustering.fit_predict(svl_matrix)
    elif args.clusterMethod == 'kmeans':
        log.debug('SVL matrix created, start clustering kmeans')

        kmeans_object = KMeans(n_clusters=args.numberOfClusters, random_state=0, n_jobs=args.threads, precompute_distances=False)
        labels_clustering = kmeans_object.fit_predict(svl_matrix)

    matrices_cluster = list(zip(matrices_list, labels_clustering))
    np.savetxt(args.outFileName, matrices_cluster, fmt="%s")
