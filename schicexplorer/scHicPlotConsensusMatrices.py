
import argparse
import os
from multiprocessing import Process, Queue
import time

import logging
log = logging.getLogger(__name__)

import cooler
from hicmatrix import HiCMatrix as hm
import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        add_help=False
    )

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--matrix', '-m',
                                help='The consensus matrix created by scHicConsensusMatrices',
                                metavar='mcool scHi-C matrix',
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the resulting cluster profile.',
                                required=False,
                                default='consensus_matrices.png')
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
    parserRequired.add_argument('--chromosomes',
                                help='List of to be plotted chromosomes',
                                nargs='+')
    return parser


def main(args=None):

    args = parse_arguments().parse_args(args)
    matrices_list = cooler.fileops.list_coolers(args.matrix)
    rows = int(np.ceil(len(matrices_list) / 4))
    columns = 4
    f, axes = plt.subplots(rows, columns)
    for i, matrix in enumerate(matrices_list):
        if args.chromosomes is not None and len(args.chromosomes) == 1:
            hic_ma = hm.hiCMatrix(pMatrixFile=args.matrix + '::' + matrix, pChrnameList=args.chromosomes)
        else:
            hic_ma = hm.hiCMatrix(pMatrixFile=args.matrix + '::' + matrix)
            if args.chromosomes:
                hic_ma.keepOnlyTheseChr(args.chromosomes)
        matrix_data = hic_ma.matrix
        matrix_data = matrix_data.toarray()
        mask = matrix_data == 0
        try:
            matrix_data[mask] = np.nanmin(matrix_data[mask == False])
        except ValueError:
            log.info('Matrix contains only 0. Set all values to {}'.format(np.finfo(float).tiny))
            matrix_data[mask] = np.finfo(float).tiny
        if np.isnan(matrix_data).any() or np.isinf(matrix_data).any():
            mask_nan = np.isnan(matrix_data)
            mask_inf = np.isinf(matrix_data)
            matrix_data[mask_nan] = np.nanmin(matrix_data[mask_nan == False])
            matrix_data[mask_inf] = np.nanmin(matrix_data[mask_inf == False])
        matrix_data += 1
        axes[i // 4, i % 4].imshow(matrix_data, cmap='RdYlBu_r', norm=LogNorm())
        axes[i // 4, i % 4].get_xaxis().set_ticks([])
        axes[i // 4, i % 4].get_yaxis().set_ticks([])

        axes[i // 4, i % 4].yaxis.set_visible(False)
        axes[i // 4, i % 4].set_xlabel(str(matrix.split('/')[-1].split('cluster_')[-1]))
    plt.savefig(args.outFileName, dpi=args.dpi)
