
import argparse
import os
from multiprocessing import Process, Queue
import time

import logging
log = logging.getLogger(__name__)

import cooler


from hicexplorer import hicAdjustMatrix
from schicexplorer._version import __version__

def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        add_help=False
    )

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='The matrix to adjust in the mcool format.',
                                required=True)
    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the adjusted matrix.',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')
    parserMutuallyExclusive = parser.add_mutually_exclusive_group()
    parserMutuallyExclusive.add_argument('--chromosomes', '-c',
                                         nargs='+',
                                         help='List of chromosomes to keep / remove')
    parserMutuallyExclusive.add_argument('--regions', '-r',
                                         help='BED file which stores a list of regions to keep / remove')
    parserMutuallyExclusive.add_argument('--maskBadRegions', '-mbr',
                                         help='Bad regions are identified and masked.')
    parserOpt.add_argument('--action',
                           help='Keep, remove or mask the list of specified chromosomes / regions ',
                           default='keep',
                           choices=['keep', 'remove', 'mask']
                           )
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser

def main(args=None):

    args = parse_arguments().parse_args(args)
    matrices_list = cooler.fileops.list_coolers(args.matrix)
    
    for i, matrix in enumerate(matrices_list):
        appendMatrix = '--appendMatrix'
        if i == 0:
            appendMatrix = ''
        args_adjust = "--matrix {} --outFileName {} --action {} {}".format(
            args.matrix + '::' + matrix,
            args.outFileName + '::' + matrix,
            args.action,
            appendMatrix
        )

        if args.chromosomes is not None:
            chromosomes = ' '.join(args.chromosomes)
            args_adjust += ' --chromosomes {}'.format(chromosomes)
        elif args.regions is not None:
            args_adjust += ' --regions {}'.format(args.regions)
        elif args.maskBadRegions is not None:
            args_adjust += ' --maskBadRegions {}'.format(args.maskBadRegions)

        args_adjust = args_adjust.split()
        log.debug('args_adjust {}'.format(args_adjust))
        hicAdjustMatrix.main(args_adjust)
