
from hicmatrix.lib import MatrixFileHandler
import argparse
from multiprocessing import Process, Queue
import time
import logging
log = logging.getLogger(__name__)
from schicexplorer._version import __version__


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Creates out of n cool files one scool file.',
        add_help=False)

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrices', '-m',
                                help='input file(s).',
                                nargs='+',
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the exported matrix.',
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

def load_cool_files(pMatricesList, pCutIntervals, pQueue):

    matrixFileHandlerList = []
    for i, matrix in enumerate(pMatricesList):
       
        matrixFileHandlerInput = MatrixFileHandler(pFileType='cool', pMatrixFile=matrix, pNoCutIntervals=True)

        _matrix, cut_intervals, nan_bins, \
            distance_counts, correction_factors = matrixFileHandlerInput.load()

        matrixFileHandlerOutput = MatrixFileHandler(pFileType='cool', pMatrixFile=matrix)

        matrixFileHandlerOutput.set_matrix_variables(_matrix,
                                                     pCutIntervals,
                                                     nan_bins,
                                                     correction_factors,
                                                     distance_counts)
        cut_intervals = None

        matrixFileHandlerList.append(matrixFileHandlerOutput)

    pQueue.put(matrixFileHandlerList)

def main(args=None):
    args = parse_arguments().parse_args(args)
    log.debug(args)
    matrix_file_handler_object_list = []
    matrixFileHandlerInput = MatrixFileHandler(pFileType='cool', pMatrixFile=args.matrices[0])

    _matrix, cut_intervals_all, nan_bins, \
        distance_counts, correction_factors = matrixFileHandlerInput.load()

  


    matrices_list = args.matrices
    # matrices_list = cell_name_list(matrices_name)
    # log.debug('size of matrices_list: {}'.format(len(matrices_list)))

    threads = args.threads

    
    matrixFileHandler_list = [None] * args.threads
    process = [None] * args.threads
    queue = [None] * args.threads

    # all_threads_done = False
    thread_done = [False] * args.threads
    # count_output = 0
    # count_call_of_read_input = 0
    # computed_pairs = 0
    matricesPerThread = len(matrices_list) // threads

    for i in range(args.threads):
        if i < threads - 1:
            matrices_name_list = matrices_list[i * matricesPerThread:(i + 1) * matricesPerThread]
        else:
            matrices_name_list = matrices_list[i * matricesPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=load_cool_files, kwargs=dict(
            # pMatrixName=matrices_name,
            pMatricesList=matrices_name_list,
            pCutIntervals=cut_intervals_all,
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
        # path_name = ''.join(matrix.split('/')[-1].split('.')[:-1])
        # matrixFileHandlerOutput.save(args.outFileName + '::/' + path_name, pApplyCorrection=True, pSymmetric=True)
