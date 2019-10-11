
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


from hicmatrix import HiCMatrix as hm
from schicexplorer.utilities import opener
from hicmatrix.lib import MatrixFileHandler

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.decomposition import TruncatedSVD

import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import LogNorm
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

    # matrices_name = args.matrix
        # threads = args.threads
    matrices_list = cooler.fileops.list_coolers(args.matrix)
    rows = int(np.ceil(len(matrices_list) / 4))
    columns =  4
    f,axes = plt.subplots(rows, columns)
    log.debug('len(axes) {}'.format(len(axes[0])))
    # if len(args.chromosomes):
        # args.chromosomes = None
    for i, matrix in enumerate(matrices_list):
        if args.chromosomes is not None and len(args.chromosomes) == 1:
            hic_ma = hm.hiCMatrix(pMatrixFile=args.matrix+ '::' +matrix, pChrnameList=args.chromosomes)
        else:
            hic_ma = hm.hiCMatrix(pMatrixFile=args.matrix+ '::' +matrix)
            if args.chromosomes:
                hic_ma.keepOnlyTheseChr(args.chromosomes)
        # matrixFileHandlerInput = MatrixFileHandler(pFileType='cool', pMatrixFile=args.matrix+ '::' +matrix, pChrnameList=args.chromosomes)
        # matrix_data, _, _, _, _ = matrixFileHandlerInput.load()
        matrix_data = hic_ma.matrix
        matrix_data = matrix_data.toarray()
        mask = matrix_data == 0
        try:
            matrix_data[mask] = np.nanmin(matrix_data[mask == False])
        except ValueError:
            log.info('Matrix contains only 0. Set all values to {}'.format(np.finfo(float).tiny))
            matrix_data[mask] = np.finfo(float).tiny
        if np.isnan(matrix_data).any() or np.isinf(matrix_data).any():
            # log.debug("any nan {}".format(np.isnan(matrix).any()))
            # log.debug("any inf {}".format(np.isinf(matrix).any()))
            mask_nan = np.isnan(matrix_data)
            mask_inf = np.isinf(matrix_data)
            matrix_data[mask_nan] = np.nanmin(matrix_data[mask_nan == False])
            matrix_data[mask_inf] = np.nanmin(matrix_data[mask_inf == False])
        matrix_data += 1
        # matrix_data = np.log(matrix_data)
        # cs = ax.imshow(z,norm=LogNorm())
        cs = axes[i//4, i%4].imshow(matrix_data, cmap='RdYlBu_r', norm=LogNorm())

        # cbar = f.colorbar(cs)

        # cbar.ax.minorticks_off()

        # # log.debug('np.array(cluster[0]).T {}'.format(np.array(cluster[0])))
        # # axes[i//4,i%4].set_yscale('log')
        # # axes[i//4, i%4].invert_yaxis()
        axes[i//4,  i%4].get_xaxis().set_ticks([])
        axes[i//4,  i%4].get_yaxis().set_ticks([])

        # if i > 0:
        axes[i//4,  i%4].yaxis.set_visible(False)
        axes[i//4,  i%4].set_xlabel(str(matrix.split('/')[-1].split('cluster_')[-1]))
    plt.savefig(args.outFileName, dpi=args.dpi)