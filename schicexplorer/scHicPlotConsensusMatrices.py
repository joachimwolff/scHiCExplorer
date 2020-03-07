
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
                                help='The consensus matrix created by scHicConsensusMatrices',
                                metavar='scool scHi-C matrix',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--outFileName', '-o',
                           help='File name to save the resulting cluster profile.',
                           required=False,
                           default='consensus_matrices.png')

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
    parserOpt.add_argument('--chromosomes', '-c',
                           help='List of to be plotted chromosomes',
                           nargs='+')
    parserOpt.add_argument('--colorMap',
                           help='Color map to use for the heatmap. Available '
                           'values can be seen here: '
                           'http://matplotlib.org/examples/color/colormaps_reference.html',
                           default='RdYlBu_r')
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def main(args=None):

    args = parse_arguments().parse_args(args)
    matrices_list = cooler.fileops.list_coolers(args.matrix)
    columns = 4
    if len(matrices_list) < columns:
        columns = len(matrices_list)
    rows = int(np.ceil(len(matrices_list) / columns))
    if rows < 1:
        rows = 1
    f, axes = plt.subplots(rows, columns)

    title_string = 'Consensus matrices of {}'.format(os.path.basename(args.matrix))
    if args.chromosomes:
        title_string += ' on chromosome: {}'.format(' '.join(args.chromosomes))
    else:
        title_string += ' on all chromosomes'
    plt.suptitle(title_string, fontsize=10)

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

        if rows == 1:
            im = axes[i % columns].imshow(matrix_data, cmap=args.colorMap, norm=LogNorm())
            axes[i % columns].get_xaxis().set_ticks([])
            axes[i % columns].get_yaxis().set_ticks([])

            axes[i % columns].yaxis.set_visible(False)
            axes[i % columns].set_xlabel(str(matrix.split('/')[-1].split('cluster_')[-1]))
        else:
            im = axes[i // columns, i % columns].imshow(matrix_data, cmap=args.colorMap, norm=LogNorm())
            axes[i // columns, i % columns].get_xaxis().set_ticks([])
            axes[i // columns, i % columns].get_yaxis().set_ticks([])

            axes[i // columns, i % columns].yaxis.set_visible(False)
            axes[i // columns, i % columns].set_xlabel(str(matrix.split('/')[-1].split('cluster_')[-1]))

    number_of_plots = len(matrices_list)
    i = -1
    while rows * columns > number_of_plots:

        axes[-1, i].axis('off')
        number_of_plots += 1
        i -= 1

    f.colorbar(im)
    plt.tight_layout()
    plt.savefig(args.outFileName, dpi=args.dpi)
    plt.close()
