
from hicmatrix.lib import MatrixFileHandler
import argparse
from multiprocessing import Process, Queue
import time
import os
import logging
log = logging.getLogger(__name__)
from schicexplorer._version import __version__
from schicexplorer.utilities import cell_name_list


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='This script offers methods to extract cool files or to update a scool file for scHiCExplorer version 5',
        add_help=False)

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='The scool matrix ',
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the exported matrix, in case of extract the folder name',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--update', '-u',
                           help='Update the scool file from the old format of scHiCExplorer until version 4 to the one used since version 5.',
                           required=False,
                           action='store_true')
    parserOpt.add_argument('--extract', '-e',
                           help='Path to a file with the cell names to be extracted into an individual cool file',
                           required=False,
                           type=str)
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module.',
                           required=False,
                           default=4,
                           type=int)
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def load_cool_files(pMatrixName, pMatricesList, pCutIntervals, pQueue):

    matrixFileHandlerList = []
    try:
        for i, matrix in enumerate(pMatricesList):

            matrixFileHandlerInput = MatrixFileHandler(pFileType='cool', pMatrixFile=pMatrixName + "::" + matrix, pNoCutIntervals=True)

            _matrix, cut_intervals, nan_bins, \
                distance_counts, correction_factors = matrixFileHandlerInput.load()

            matrixFileHandlerOutput = MatrixFileHandler(pFileType='cool', pMatrixFile=matrix.split('/')[-1])

            matrixFileHandlerOutput.set_matrix_variables(_matrix,
                                                         pCutIntervals,
                                                         nan_bins,
                                                         correction_factors,
                                                         distance_counts)
            cut_intervals = None

            matrixFileHandlerList.append(matrixFileHandlerOutput)
    except Exception as exp:
        pQueue.put('Fail: ' + str(exp) + traceback.format_exc())
        return
    pQueue.put(matrixFileHandlerList)


def main(args=None):
    args = parse_arguments().parse_args(args)
    log.debug(args)
    matrix_file_handler_object_list = []

    matrices_list = cell_name_list(args.matrix)
    if args.extract is not None:
        matrix_list_tmp = []
        with open(args.extract, 'r') as file:
            for line in file:
                values = line.strip()
                if not values.startswith('/cells'):
                    values = '/cells/' + values
                if values in matrices_list:
                    matrix_list_tmp.append(values)

        matrices_list = matrix_list_tmp

    if len(matrices_list) == 0:
        raise OSError('No cells for processing. Terminating.')
        exit(1)
    if len(matrices_list) < args.threads:
        args.threads = len(matrices_list)

    matrixFileHandlerInput = MatrixFileHandler(pFileType='cool', pMatrixFile=args.matrix + "::" + matrices_list[0])

    _matrix, cut_intervals_all, nan_bins, \
        distance_counts, correction_factors = matrixFileHandlerInput.load()

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
        process[i] = Process(target=load_cool_files, kwargs=dict(
            pMatrixName=args.matrix,
            pMatricesList=matrices_name_list,
            pCutIntervals=cut_intervals_all,
            pQueue=queue[i]
        )
        )
        process[i].start()

    all_data_collected = False
    fail_flag = False
    fail_message = ''
    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                matrixFileHandler_list[i] = queue[i].get()
                if 'Fail:' in matrixFileHandler_list[i]:
                    fail_flag = True
                    fail_message = matrixFileHandler_list[i][6:]
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

    if fail_flag:
        log.error(fail_message)
        exit(1)
    matrix_file_handler_object_list = [item for sublist in matrixFileHandler_list for item in sublist]

    if args.update:
        matrixFileHandler = MatrixFileHandler(pFileType='scool')
        matrixFileHandler.matrixFile.coolObjectsList = matrix_file_handler_object_list
        matrixFileHandler.save(args.outFileName, pSymmetric=True, pApplyCorrection=False)
    else:
        if not os.path.exists(args.outFileName):
            try:
                os.makedirs(args.outFileName)
            except OSError as exc:  # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise
        for matrixFileHandler in matrix_file_handler_object_list:
            matrixFileHandler.save(args.outFileName + '/' + matrixFileHandler.matrixFile.matrixFileName + '.cool', pApplyCorrection=True, pSymmetric=True)
