import argparse
from multiprocessing import Process, Queue
import time
import logging
log = logging.getLogger(__name__)

import cooler
import numpy as np

from hicmatrix import HiCMatrix
from hicmatrix.lib import MatrixFileHandler

def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        add_help=False
    )

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--matrix', '-m',
                                help='The single cell Hi-C interaction matrices to investigate for QC. Needs to be in mcool format',
                                metavar='mcool scHi-C matrix',
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name of the consensus mcool matrix.',
                                required=True)

    parserRequired.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module.',
                           required=False,
                           default=4,
                           type=int)

    return parser

def compute_sum(pMatrixName, pMatricesList, pThread, pQueue):
    sum_list = []
    for i, matrix in enumerate(pMatricesList):
       
        matrixFileHandler = MatrixFileHandler(pFileType='cool', pMatrixFile=pMatrixName+ '::' +matrix)
        _matrix, cut_intervals, nan_bins, \
                distance_counts, correction_factors = matrixFileHandler.load()
        try:
            sum_of_matrix = _matrix.sum()
        except:
            sum_list.append()
        sum_list.append(sum_of_matrix)
    pQueue.put(sum_list)

def compute_normalize(pMatrixName, pMatricesList, pArgminSum, pSumOfAll, pAppend, pQueue):
   
    matrixFileHandlerList = []
    for i, matrix in enumerate(pMatricesList):
        if i == 0 and pAppend:
            append = True
        else:
            append = False
        matrixFileHandler = MatrixFileHandler(pFileType='cool', pMatrixFile=pMatrixName+ '::' +matrix)
        _matrix, cut_intervals, nan_bins, \
                distance_counts, correction_factors = matrixFileHandler.load()
        _matrix.data = _matrix.data.astype(np.float32)
        # if i != pArgmin:
        mask = np.isnan(_matrix.data)
        _matrix.data[mask] = 0

        mask = np.isinf(_matrix.data)
        _matrix.data[mask] = 0
        adjust_factor = pSumOfAll[i] / pArgminSum
        _matrix.data /= adjust_factor
        mask = np.isnan(_matrix.data)

        mask = np.isnan(_matrix.data)
        _matrix.data[mask] = 0

        mask = np.isinf(_matrix.data)
        _matrix.data[mask] = 0
        _matrix.eliminate_zeros()

        matrixFileHandlerOutput = MatrixFileHandler(pFileType='cool', pAppend=pAppend, pEnforceInteger=False, pFileWasH5=False, pHic2CoolVersion=None)

        matrixFileHandlerOutput.set_matrix_variables(_matrix, cut_intervals, nan_bins,
                                                    correction_factors, distance_counts)
        # matrixFileHandler.set_matrix_variables(_matrix, cut_intervals, nan_bins,
        #                                             correction_factors, distance_counts)
    
        matrixFileHandlerList.append(matrixFileHandlerOutput)
    
    pQueue.put(matrixFileHandlerList)

def main(args=None):

    args = parse_arguments().parse_args(args)

    matrices_name = args.matrix
    matrices_list = cooler.fileops.list_coolers(matrices_name)
    
    threads = args.threads

    sum_list_threads = [None] * args.threads
    process = [None] * args.threads
    all_data_processed = False
    queue = [None] * args.threads

    all_threads_done = False
    thread_done = [False] * args.threads
    count_output = 0
    count_call_of_read_input = 0
    computed_pairs = 0
    matricesPerThread = len(matrices_list) // threads

    for i in range(args.threads):
        if i < threads - 1:
            matrices_name_list = matrices_list[i * matricesPerThread:(i + 1) * matricesPerThread]
        else:
            matrices_name_list = matrices_list[i * matricesPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=compute_sum, kwargs=dict(
                            pMatrixName=matrices_name,
                            pMatricesList=matrices_name_list, 
                            pThread=i,
                            pQueue=queue[i]
            )
        )
        process[i].start()
                    
    all_data_collected = False
    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                sum_list_threads[i] = queue[i].get()
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

    sum_of_all = [item for sublist in sum_list_threads for item in sublist]
    
    argmin = np.argmin(sum_of_all)
    argminSum = sum_of_all[argmin]

    matricesPerThread = len(matrices_list) // threads

    matrixFileHandlerListThreads = [None] * args.threads
    process = [None] * args.threads
    all_data_processed = False
    queue = [None] * args.threads

    all_threads_done = False
    thread_done = [False] * args.threads
    count_output = 0
    count_call_of_read_input = 0
    computed_pairs = 0

    for i in range(args.threads):
        if i < args.threads - 1:
            matrices_name_list = matrices_list[i * matricesPerThread:(i + 1) * matricesPerThread]
            sum_of_all_list = sum_of_all[i * matricesPerThread:(i + 1) * matricesPerThread]
        else:
            matrices_name_list = matrices_list[i * matricesPerThread:]
            sum_of_all_list = sum_of_all[i * matricesPerThread:]

        log.debug('thread: {} size of matrix list {}, size of sum list {}'.format(i, len(matrices_name_list), len(sum_of_all_list)))
        queue[i] = Queue()
        process[i] = Process(target=compute_normalize, kwargs=dict(
                            pMatrixName = matrices_name,
                            pMatricesList= matrices_name_list, 
                            pArgminSum=argminSum,
                            pSumOfAll=sum_of_all_list,
                            pAppend = i > 0,
                            pQueue=queue[i]
            )
        )
        process[i].start()
    all_data_collected = False
    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                matrixFileHandlerListThreads[i] = queue[i].get()

                queue[i] = None
                process[i].join()
                process[i].terminate()
                process[i] = None
                thread_done[i] = True
                log.debug('Data collected: {}'.format(i))
                log.debug('Data collected: {}'.format(thread_done))

        all_data_collected = True
        for thread in thread_done:
            if not thread:
                all_data_collected = False
        time.sleep(1)
  
    matrixFileHandlerList = [item for sublist in matrixFileHandlerListThreads for item in sublist]


    for i, matrixFileHandler in enumerate(matrixFileHandlerList):
        matrixFileHandler.save(args.outFileName + '::' + matrices_list[i] , pSymmetric=True, pApplyCorrection=False)
        
