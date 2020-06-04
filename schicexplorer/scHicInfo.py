import argparse
from multiprocessing import Process, Queue
import time
import logging
log = logging.getLogger(__name__)

import cooler
import numpy as np

from hicmatrix.lib import MatrixFileHandler
from schicexplorer._version import __version__
from schicexplorer.utilities import cell_name_list



def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description=''
    )

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='The single cell Hi-C interaction matrices to investigate for QC. Needs to be in scool format',
                                metavar='scool scHi-C matrix',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def main(args=None):

    args = parse_arguments().parse_args(args)

    matrices_list = cell_name_list(args.matrix)

    print('Filename: {}'.format(args.matrix))
    print('Contains {} single-cell matrices'.format(len(matrices_list)))
    print('The information stored via cooler.info of the first cell is: \n')
    cooler_file = cooler.Cooler(args.matrix + '::' + matrices_list[0])

    if cooler_file.info is not None:
        for key, value in cooler_file.info.items():
            print(key, value)
