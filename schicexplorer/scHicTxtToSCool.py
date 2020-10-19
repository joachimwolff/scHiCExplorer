
from hicmatrix.lib import MatrixFileHandler
import argparse
from multiprocessing import Process, Queue
import time
import logging

from scipy.sparse import csr_matrix
import numpy as np

log = logging.getLogger(__name__)
from schicexplorer._version import __version__


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Creates out of n txt files from Ramani 2017 one scool file. ',
        add_help=False)

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrices', '-m',
                                help='input file(s).',
                                nargs='+',
                                required=True)
    parserRequired.add_argument('--chromosomeSizes', '-cs',
                                help='Text file with two columns, first column is the name of the chromosome, second one its length.',
                                required=True)
    parserRequired.add_argument('--resolution', '-r',
                                help='The resolution of the matrix.',
                                required=True,
                                type=int)
    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the scool matrix.',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module.',
                           required=False,
                           default=4,
                           type=int)
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def txt_to_matrixFileHandler(pMatricesList, pMatrixDimensions, pCutIntervals, pQueue):

    matrixFileHandlerList = []

    # matrix_dimensions = pGenomeLength // pResolution
    for i, matrix in enumerate(pMatricesList):

        
    
        # create csr matrix
        # hic_matrix = csr_matrix((pMatrixDimensions,pMatrixDimensions), dtype=float)
        instances = []
        features = []
        data = []
        with open(matrix, 'r') as file:
            for i, line in enumerate(file.readlines()):
                line = line.strip()
                if len(line) == 0:
                    continue
                x, y, count = line.split('\t')[:3]
                instances.append(int(x))
                features.append(int(y))
                data.append(float(count))

        cell_type = matrix.split('_')[2]

        # nan_bins, correction-factors, dostance_counts = None
        log.debug('matrix name {}'.format(matrix))

        log.debug('max(instances) {} max(features) {} pMatrixDimensions {}'.format(max(instances),max(features), pMatrixDimensions))
        hic_matrix = csr_matrix((data, (instances, features)), (pMatrixDimensions, pMatrixDimensions), dtype=np.float)

        # matrixFileHandlerInput = MatrixFileHandler(pFileType='cool', pMatrixFile=matrix)

        # _matrix, cut_intervals, nan_bins, \
        #     distance_counts, correction_factors = matrixFileHandlerInput.load()

        matrixFileHandlerOutput = MatrixFileHandler(pFileType='cool', pMatrixFile=matrix)

        matrixFileHandlerOutput.set_matrix_variables(hic_matrix,
                                                     pCutIntervals,
                                                     None,
                                                     None,
                                                     None)
                                                
        if matrixFileHandlerOutput.matrixFile.hic_metadata is None:
            matrixFileHandlerOutput.matrixFile.hic_metadata = {}
            matrixFileHandlerOutput.matrixFile.hic_metadata['cell_type'] = cell_type

        matrixFileHandlerList.append(matrixFileHandlerOutput)

    pQueue.put(matrixFileHandlerList)


def main(args=None):
    args = parse_arguments().parse_args(args)
    log.debug(args)
    matrix_file_handler_object_list = []
    # matrixFileHandlerInput = MatrixFileHandler(pFileType='cool', pMatrixFile=args.matrices[0])

    # _matrix, cut_intervals_all, nan_bins, \
    #     distance_counts, correction_factors = matrixFileHandlerInput.load()


    # read genome sizes
    chromosome_sizes = {}
    genome_size = 0
    matrix_dimensions = 0
    with open(args.chromosomeSizes, 'r') as file:
        for i, line in enumerate(file.readlines()):
            line = line.strip()
            chromosome_name, chromosome_size = line.split('\t')
            chromosome_sizes[chromosome_name] = int(chromosome_size)
            genome_size += int(chromosome_size)
            matrix_dimensions += ((int(chromosome_size) // args.resolution) + 1)
    
    log.debug('genome_size {}'.format(genome_size))
    log.debug('args.resolution {}'.format(args.resolution))


    # matrix_dimensions = genome_size // args.resolution
    log.debug('matrix_dimensions {}'.format(matrix_dimensions))

    log.debug('chromosome_sizes {}'.format(chromosome_sizes))
    # create cut_intervals:
    cut_intervals = []
    
    for chromosome, size in chromosome_sizes.items():
        for interval in range(0, size, args.resolution):
            cut_intervals.append((chromosome, interval,
                                min(size, interval + args.resolution), 1))

    matrices_list = args.matrices
    
    threads = args.threads

    matrixFileHandler_list = [None] * args.threads
    process = [None] * args.threads
    queue = [None] * args.threads

    thread_done = [False] * args.threads
    matricesPerThread = len(matrices_list) // threads

    for i in range(args.threads):
        if i < threads - 1:
            matrices_name_list = matrices_list[i * matricesPerThread:(i + 1) * matricesPerThread]
        else:
            matrices_name_list = matrices_list[i * matricesPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=txt_to_matrixFileHandler, kwargs=dict(
            pMatricesList=matrices_name_list,
            pMatrixDimensions=matrix_dimensions,
            pCutIntervals=cut_intervals,
            pQueue=queue[i]
        )
        )
        process[i].start()

    all_data_collected = False
    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                matrixFileHandler_list[i] = queue[i].get()
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

    matrix_file_handler_object_list = [item for sublist in matrixFileHandler_list for item in sublist]

    matrixFileHandler = MatrixFileHandler(pFileType='scool')
    matrixFileHandler.matrixFile.coolObjectsList = matrix_file_handler_object_list
    matrixFileHandler.save(args.outFileName, pSymmetric=True, pApplyCorrection=False)
