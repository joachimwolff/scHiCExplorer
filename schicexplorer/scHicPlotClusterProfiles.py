# read all matrices
# get non-zeros and flatten it (x*length) + y
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
logging.getLogger('hicmatrix').setLevel(logging.ERROR)


log = logging.getLogger(__name__)

import cooler

from hicmatrix import HiCMatrix as hm
from hicmatrix.lib import MatrixFileHandler

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from multiprocessing import Process, Queue


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
    read_coverage = []
    sparsity = []
    read_distribution = []
    # binSize = 1000000
    length_old = 50
    for i, matrix in enumerate(pMatricesList):
        if pChromosomes is not None and len(pChromosomes) == 1:
            hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName + '::' + matrix, pChrnameList=pChromosomes)
        else:
            hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName + '::' + matrix)
            if pChromosomes:
                hic_ma.keepOnlyTheseChr(pChromosomes)
        _matrix = hic_ma.matrix
        maximalDistance = pMaximalDistance // hic_ma.getBinSize()

        instances, features = _matrix.nonzero()

        relative_distance = np.absolute(instances - features)
        read_distribution_ = np.zeros(maximalDistance)
        for j, relative_distance_ in enumerate(relative_distance):
            if relative_distance_ < maximalDistance:
                read_distribution_[relative_distance_] += _matrix.data[j]
        read_distribution_ /= np.sum(_matrix.data)
        read_distribution.append(read_distribution_)

    pQueue.put(read_distribution)


def main(args=None):

    args = parse_arguments().parse_args(args)

    # if args.loadFromNpz:

    matrices_name = args.matrix
    threads = args.threads
    matrices_list = cooler.fileops.list_coolers(matrices_name)
    read_coverage = [None] * threads
    sparsity = [None] * threads

    all_data_collected = False
    thread_done = [False] * threads
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
                read_coverage[i] = queue[i].get()

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
    length_old = len(read_coverage[0][0])
    for thread_data in read_coverage:
        for matrix_data in thread_data:
            read_distributions.append(matrix_data)

    read_distributions = np.array(read_distributions)

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

    # w=10
    # h=100
    columns = len(clusters)
    rows = 1
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

    fig, axes = plt.subplots(rows, columns,
                             gridspec_kw={'width_ratios': cluster_size})
    for i, cluster in enumerate(clusters_list):
        try:
            im = axes[i].imshow(cluster.T, cmap='RdYlBu_r', norm=LogNorm(), aspect="auto")
            axes[i].invert_yaxis()
            axes[i].get_xaxis().set_ticks([])
            if i > 0:
                axes[i].yaxis.set_visible(False)
            axes[i].set_title('{}'.format(len(cluster)), fontsize=8)
            axes[i].set_xlabel(str(i))
        except Exception:
            log.debug('Cluster {} failed. Elements: {}'.format(i, len(cluster)))
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    plt.savefig(args.outFileName, dpi=args.dpi)

    plt.close()
