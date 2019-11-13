import argparse
from multiprocessing import Process, Queue
import time
import logging
log = logging.getLogger(__name__)

import cooler
import numpy as np

from hicmatrix.lib import MatrixFileHandler


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        add_help=False
    )

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='The single cell Hi-C interaction matrices to investigate for QC. Needs to be in mcool format',
                                metavar='mcool scHi-C matrix',
                                required=True)

    return parser


def main(args=None):

    args = parse_arguments().parse_args(args)

    matrices_list = cooler.fileops.list_coolers(args.matrix)

    print('Filename: {}'.format(args.matrix))
    print('Contains {} single-cell matrices'.format(len(matrices_list)))

    
    
