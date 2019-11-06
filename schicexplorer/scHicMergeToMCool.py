
from hicmatrix.lib import MatrixFileHandler
import argparse

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Creates out of n cool files one mcool file.',
        add_help=False)

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrices', '-m',
                                help='input file(s).',
                                nargs='+',
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the exported matrix.',
                                required=True)
    return parser


def main(args=None):
    args = parse_arguments().parse_args(args)
    log.debug(args)

    for i, matrix in enumerate(args.matrices):
        matrixFileHandlerInput = MatrixFileHandler(pFileType='cool', pMatrixFile=matrix)

        _matrix, cut_intervals, nan_bins, \
            distance_counts, correction_factors = matrixFileHandlerInput.load()

        append = False
        if i > 0:
            append = True
        matrixFileHandlerOutput = MatrixFileHandler(pFileType='cool', pAppend=append)

        matrixFileHandlerOutput.set_matrix_variables(_matrix,
                                                     cut_intervals,
                                                     nan_bins,
                                                     correction_factors,
                                                     distance_counts)

        path_name = ''.join(matrix.split('/')[-1].split('.')[:-1])
        matrixFileHandlerOutput.save(args.outFileName + '::/' + path_name, pApplyCorrection=True, pSymmetric=True)
