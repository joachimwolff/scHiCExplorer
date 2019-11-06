import argparse
import os
from multiprocessing import Process, Queue
import time

import logging
log = logging.getLogger(__name__)

import cooler

from hicmatrix import HiCMatrix as hm
from hicmatrix.lib import MatrixFileHandler

import numpy as np


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


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
                                help='Text file which contains per matrix the associated cluster.',
                                metavar='cluster file',
                                required=True)
    parserRequired.add_argument('--chromosomes',
                                help='List of to be plotted chromosomes',
                                nargs='+')
    parserRequired.add_argument('--maximalDistance', '-md',
                                help='maximal distance in bases',
                                required=False,
                                default=50000000,
                                type=int)
    parserRequired.add_argument('--shortRange', '-sr',
                                help='Distance to split short and long range contacts to sort the samples per cluster.',
                                required=False,
                                default=20000000,
                                type=int)
    parserRequired.add_argument('--orderBy', '-ob',
                                help='Algorithm to cluster the Hi-C matrices',
                                choices=['svl', 'orderByFile'],
                                default='spectral')
    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the resulting cluster profile.',
                                required=False,
                                default='clusters_profile.png')
    parserRequired.add_argument('--dpi', '-d',
                                help='The dpi of the plot.',
                                required=False,
                                default=300,
                                type=int)
    parserRequired.add_argument('--threads',
                                help='Number of threads. Using the python multiprocessing module.',
                                required=False,
                                default=4,
                                type=int)
    return parser


def compute_read_distribution(pMatrixName, pMatricesList, pMaximalDistance, pChromosomes, pQueue):
    read_distribution = []
    resolution = 0
    for i, matrix in enumerate(pMatricesList):
        if pChromosomes is not None and len(pChromosomes) == 1:
            hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName + '::' + matrix, pChrnameList=pChromosomes)
        else:
            hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName + '::' + matrix)
            if pChromosomes:
                hic_ma.keepOnlyTheseChr(pChromosomes)
        _matrix = hic_ma.matrix
        resolution = hic_ma.getBinSize()
        maximalDistance = pMaximalDistance // resolution

        instances, features = _matrix.nonzero()

        relative_distance = np.absolute(instances - features)
        read_distribution_ = np.zeros(maximalDistance)
        sum_of_matrix_within_maximalDistance = 0
        for j, relative_distance_ in enumerate(relative_distance):
            if relative_distance_ < maximalDistance:
                read_distribution_[relative_distance_] += _matrix.data[j]
                sum_of_matrix_within_maximalDistance += _matrix.data[j]
        read_distribution_ /= sum_of_matrix_within_maximalDistance
        read_distribution.append(read_distribution_)

    pQueue.put([read_distribution, resolution])


def main(args=None):

    args = parse_arguments().parse_args(args)

    matrices_name = args.matrix
    threads = args.threads
    matrices_list = cooler.fileops.list_coolers(matrices_name)
    read_coverage = [None] * threads

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
        process[i] = Process(target=compute_read_distribution, kwargs=dict(
            pMatrixName=matrices_name,
            pMatricesList=matrices_name_list,
            pMaximalDistance=args.maximalDistance,
            pChromosomes=args.chromosomes,
            pQueue=queue[i]
        )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                read_coverage[i], resolution = queue[i].get()

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

    read_distributions = []
    for thread_data in read_coverage:
        for matrix_data in thread_data:
            read_distributions.append(matrix_data)

    read_distributions = np.array(read_distributions)

    clusters = {}
    clusters_svl = {}
    # resolution =
    short_range_distance = args.shortRange // resolution
    mitotic_range_distance = 12000000 // resolution

    with open(args.clusters, 'r') as cluster_file:

        for i, line in enumerate(cluster_file.readlines()):
            line = line.strip()
            line_ = line.split(' ')[1]
            if int(line_) in clusters:
                clusters[int(line_)].append(read_distributions[i])
                clusters_svl[int(line_)].append(np.sum(read_distributions[i][:short_range_distance]) / np.sum(read_distributions[i][short_range_distance:mitotic_range_distance]))
            else:
                clusters[int(line_)] = [read_distributions[i]]
                clusters_svl[int(line_)] = [np.sum(read_distributions[i][:short_range_distance]) / np.sum(read_distributions[i][short_range_distance:mitotic_range_distance])]
    if args.orderBy == 'svl':
        for i, cluster_key in enumerate(clusters.keys()):
            clusters[cluster_key] = np.array(clusters[cluster_key])
            clusters_svl[cluster_key] = np.array(clusters_svl[cluster_key])
            sorted_indices = np.argsort(clusters_svl[cluster_key])
            clusters[cluster_key] = clusters[cluster_key][sorted_indices]

    cluster_to_plot = []

    clusters_list = []
    cluster_size = []
    for i, key in enumerate(clusters):
        cluster_to_plot = []

        for cluster_item in clusters[key]:
            cluster_to_plot.append(cluster_item)
        clusters_list.append(np.array(cluster_to_plot))
        cluster_size.append(len(cluster_to_plot))

    cluster_size = np.array(cluster_size)
    cluster_size = (1.0 - 0.1) * (cluster_size - min(cluster_size)) / (max(cluster_size) - min(cluster_size)) + (0.1)
    cluster_size = list(cluster_size)

    all_data = None
    index_clusters = []
    cluster_ticks = []
    ticks_position = []
    for i, cluster in enumerate(clusters_list):
        if all_data is None:
            all_data = cluster
            index_clusters.append(len(cluster))
            ticks_position.append(0 + len(cluster) // 2)
        else:
            all_data = np.append(all_data, cluster, axis=0)
            index_clusters.append(index_clusters[i - 1] + len(cluster))
            ticks_position.append(index_clusters[i - 1] + len(cluster) // 2)

        cluster_ticks.append('Cluster {} ({} cells)'.format(i + 1, len(cluster)))

    fig = plt.figure(figsize=(10, 2))

    plt.imshow(all_data.T, cmap='RdYlBu_r', norm=LogNorm(), aspect="auto")

    for index in index_clusters[:-1]:
        plt.axvline(index, color='black', linewidth=0.4)

    y_ticks = []
    y_labels = []
    for i in range(0, args.maximalDistance // resolution, 1):
        if i % 10 == 0:
            y_ticks.append(i)

            y_labels.append(str(i) + 'MB')
    plt.yticks(ticks=y_ticks, labels=y_labels, fontsize=8)

    plt.gca().invert_yaxis()
    plt.xticks(ticks=ticks_position, labels=cluster_ticks, fontsize=4)
    fig.autofmt_xdate()

    cbar = plt.colorbar()
    cbar.ax.set_ylabel('% contacts', rotation=270, fontsize=5)
    cbar.ax.invert_yaxis()
    cbar.ax.tick_params(labelsize=5)
    plt.tight_layout()
    plt.savefig(args.outFileName, dpi=args.dpi)

    plt.close()
