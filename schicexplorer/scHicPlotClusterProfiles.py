import argparse
import os
from multiprocessing import Process, Queue
import time

import logging
log = logging.getLogger(__name__)

import cooler

from hicmatrix import HiCMatrix as hm
from hicmatrix.lib import MatrixFileHandler
from schicexplorer.utilities import cell_name_list

import numpy as np


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from schicexplorer._version import __version__


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description=''
    )

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--matrix', '-m',
                                help='The single cell Hi-C interaction matrices to investigate for QC. Needs to be in scool format',
                                metavar='scool scHi-C matrix',
                                required=True)
    parserRequired.add_argument('--clusters', '-c',
                                help='Text file which contains per matrix the associated cluster.',
                                metavar='cluster file',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--chromosomes',
                           help='List of to be plotted chromosomes',
                           nargs='+')
    parserOpt.add_argument('--maximalDistance', '-md',
                           help='maximal distance in bases',
                           required=False,
                           default=50000000,
                           type=int)
    parserOpt.add_argument('--distanceShortRange', '-ds',
                           help='Distance which should be considered as short range. Default 2MB.',
                           default=2000000,
                           type=int)
    parserOpt.add_argument('--distanceLongRange', '-dl',
                           help='Distance which should be considered as short range. Default 12MB.',
                           default=12000000,
                           type=int)
    parserOpt.add_argument('--orderBy', '-ob',
                           help='Algorithm to cluster the Hi-C matrices',
                           choices=['svl', 'orderByFile'],
                           default='svl')
    parserOpt.add_argument('--fontsize',
                           help='Fontsize in the plot for x and y axis.',
                           type=float,
                           default=10)
    parserOpt.add_argument('--rotationX',
                           help='Rotation in degrees for the labels of x axis.',
                           type=int,
                           default=0)
    parserMutuallyExclusiveGroupFilter = parser.add_mutually_exclusive_group(
        required=False)
    parserMutuallyExclusiveGroupFilter.add_argument('--ticks',
                           help='Plot the cluster names as ticks. Use legend if they overlap.',
                           action='store_true')
    parserMutuallyExclusiveGroupFilter.add_argument('--legend',
                           help='Plot the cluster names as legend.  Might be helpful if the ticks overlap.',
                           action='store_true')
    parserOpt.add_argument('--outFileName', '-o',
                           help='File name to save the resulting cluster profile.',
                           required=False,
                           default='clusters_profile.png')
    parserOpt.add_argument('--dpi', '-d',
                           help='The dpi of the plot.',
                           required=False,
                           default=300,
                           type=int)
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module.',
                           required=False,
                           default=4,
                           type=int)
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def compute_read_distribution(pMatrixName, pMatricesList, pMaximalDistance, pChromosomes, pQueue):
    read_distribution = []
    # resolution = 0
    # resolution = hic_ma.getBinSize()

    hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName + '::' + pMatricesList[0])
    resolution = hic_ma.getBinSize()

    for i, matrix in enumerate(pMatricesList):
        # if pChromosomes is not None and len(pChromosomes) == 1:
        #     hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName + '::' + matrix, pChrnameList=pChromosomes, )
        # else:
        #     hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName + '::' + matrix)
        #     if pChromosomes:
        #         hic_ma.keepOnlyTheseChr(pChromosomes)
        # _matrix = hic_ma.matrix
        # log.debug('name of matrix: {}'.format(pMatrixName + '::' + matrix))
        if pChromosomes is not None and len(pChromosomes) == 1:
            hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName + '::' + matrix, pChrnameList=pChromosomes, pNoIntervalTree=True, pUpperTriangleOnly=True, pLoadMatrixOnly=True, pRestoreMaskedBins=False)
        else:
            if not pChromosomes:
                hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName + '::' + matrix, pNoIntervalTree=True, pUpperTriangleOnly=True, pLoadMatrixOnly=True, pRestoreMaskedBins=False)
            else:
                hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName + '::' + matrix, pNoIntervalTree=False, pUpperTriangleOnly=True, pLoadMatrixOnly=True, pRestoreMaskedBins=False)
            if pChromosomes:
                hic_ma.keepOnlyTheseChr(pChromosomes)
        _matrix = hic_ma.matrix
        # resolution = hic_ma.getBinSize()
        maximalDistance = pMaximalDistance // resolution

        # instances, features = _matrix.nonzero()
        instances = _matrix[0]
        features = _matrix[1]
        data = _matrix[2]

        relative_distance = np.absolute(instances - features)
        read_distribution_ = np.zeros(maximalDistance)
        sum_of_matrix_within_maximalDistance = 0
        for j, relative_distance_ in enumerate(relative_distance):
            if relative_distance_ < maximalDistance:
                read_distribution_[relative_distance_] += data[j]
                sum_of_matrix_within_maximalDistance += data[j]
        read_distribution_ /= sum_of_matrix_within_maximalDistance
        read_distribution.append(read_distribution_)

    pQueue.put([read_distribution, resolution])


def main(args=None):

    args = parse_arguments().parse_args(args)

    matrices_name = args.matrix
    threads = args.threads
    matrices_list = cell_name_list(matrices_name)
    read_coverage = [None] * threads

    all_data_collected = False
    thread_done = [False] * threads
    length_index = [None] * threads
    length_index[0] = 0
    # matrices_list = matrices_list[:16]
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
    short_range_distance = args.distanceShortRange // resolution
    long_range_distance = args.distanceLongRange // resolution

    with open(args.clusters, 'r') as cluster_file:

        for i, line in enumerate(cluster_file.readlines()):
            line = line.strip()
            line_ = line.split(' ')[1]
            if int(line_) in clusters:
                clusters[int(line_)].append(read_distributions[i])
                clusters_svl[int(line_)].append(np.sum(read_distributions[i][:short_range_distance]) / np.sum(read_distributions[i][short_range_distance:long_range_distance]))
            else:
                clusters[int(line_)] = [read_distributions[i]]
                clusters_svl[int(line_)] = [np.sum(read_distributions[i][:short_range_distance]) / np.sum(read_distributions[i][short_range_distance:long_range_distance])]
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
    cluster_ticks_top = []
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

        cluster_ticks.append('Cluster {}: {} cells'.format((i+1), len(cluster)))
        cluster_ticks_top.append('Cluster {}'.format(i + 1))


    if len(matrices_list) > 1000:
        fig = plt.figure(figsize=(10, 3))
    elif len(matrices_list) > 5000:
        fig = plt.figure(figsize=(5, 3))
    elif len(matrices_list) > 250:
        fig = plt.figure(figsize=(4, 3))
    else:
        fig = plt.figure(figsize=(3, 3))

    plt.imshow(all_data.T, cmap='RdYlBu_r', norm=LogNorm(), aspect="auto")

    for index in index_clusters[:-1]:
        plt.axvline(index, color='black', linewidth=0.4)

    y_ticks = []
    y_labels = []

    log.debug('args.maximalDistance {} '.format(args.maximalDistance))
    unit = 'MB'
    if args.maximalDistance <= 1000:
        factor = 100
        unit = 'kB'

    elif args.maximalDistance <= 10000:
        factor = 1000
    elif args.maximalDistance <= 100000:
        factor = 10000
    elif args.maximalDistance <= 1000000:
        factor = 100000
    elif args.maximalDistance <= 10000000:
        factor = 1000000
    elif args.maximalDistance <= 100000000:
        factor = 10000000
    else:
        factor = 100000000
    # resolution = 1000000
    for i in range(0, (args.maximalDistance) + 1, resolution):
        if i % (factor) == 0:
            log.debug('i {} resolution {}'.format(i // resolution, resolution))
            y_ticks.append(i // resolution)

            y_labels.append(str(i // factor * 10) + unit)
        else:
            log.debug('i {} resolution {} factor {}'.format(i, resolution, factor))

    plt.yticks(ticks=y_ticks, labels=y_labels, fontsize=args.fontsize)

    plt.gca().invert_yaxis()
    if args.ticks:
        plt.xticks(ticks=ticks_position, labels=cluster_ticks_top, rotation=args.rotationX, fontsize=args.fontsize)

        
        # secax = plt.secondary_xaxis('top', functions=(deg2rad, rad2deg))
        # secax.set_xlabel('angle [rad]')
        # plt.xticks(ticks=ticks_position, labels=cluster_ticks_top, rotation=args.rotationX, fontsize=args.fontsize)

    elif args.legend:
        plt.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False)
        leg = plt.legend(cluster_ticks, loc='upper center', bbox_to_anchor=(0.5, -0.5),
          fancybox=True, shadow=False, ncol=3, fontsize=args.fontsize)
        for item in leg.legendHandles:
            item.set_visible(False)
    else:
        plt.tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                labelbottom=False) 
    fig.autofmt_xdate()

    cbar = plt.colorbar()
    cbar.ax.set_ylabel('% contacts', rotation=270, fontsize=args.fontsize)
    cbar.ax.invert_yaxis()
    cbar.ax.tick_params(labelsize=5)
    
    
    # plt.text(0.05, 0.95, cluster_ticks, fontsize=args.fontsize,
    #             bbox_to_anchor=(0.5, -0.5))
    # props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    # plt.figtext(0, 0, cluster_ticks, wrap=True,
    #         horizontalalignment='center', fontsize=args.fontsize, bbox=props)

    # plt.annotate(['Label', 'label2'], xy=(-12, -12), xycoords='axes points',
    #         size=14, ha='right', va='top',
    #         bbox=dict(boxstyle='round', fc='w'))
    plt.tight_layout()
    plt.savefig(args.outFileName, dpi=args.dpi)

    plt.close()
