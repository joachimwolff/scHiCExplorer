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
    parserRequired.add_argument('--clusters', '-c',
                                help='Text file which contains per matrix the assoziated cluster.',
                                metavar='cluster file',
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


def compute_consensus_matrix(pMatrixName, pMatricesList, pAppend, pCounter, pQueue):
    read_coverage = []
    sparsity = []
    read_distribution = []
    max_length = 0
    consensus_matrix = None
    for i, matrix in enumerate(pMatricesList):
        matrixFileHandlerInput = MatrixFileHandler(pFileType='cool', pMatrixFile=pMatrixName+ '::' +matrix)
        _matrix, cut_intervals, nan_bins, \
                distance_counts, correction_factors = matrixFileHandlerInput.load()
    
        if consensus_matrix is None:
            consensus_matrix = _matrix
        else:
            consensus_matrix += _matrix
    
    hic2CoolVersion = matrixFileHandlerInput.matrixFile.hic2cool_version
    matrixFileHandlerOutput = MatrixFileHandler(pFileType='cool', pAppend=pAppend, pEnforceInteger=False, pFileWasH5=False, pHic2CoolVersion=hic2CoolVersion)

    matrixFileHandlerOutput.set_matrix_variables(consensus_matrix, cut_intervals, nan_bins,
                                                    correction_factors, distance_counts)
    pQueue.put([matrixFileHandlerOutput, pCounter])

def main(args=None):

    args = parse_arguments().parse_args(args)

    matrices_name = args.matrix
    matrices_list = cooler.fileops.list_coolers(matrices_name)
    clusters = {}
    with open(args.clusters, 'r') as cluster_file:

        for i, line in enumerate(cluster_file.readlines()):
            line = line.strip()
            file_path, cluster = line.split(' ')
            
            if int(cluster) in clusters:
                clusters[int(cluster)].append(file_path)
            else:
                clusters[int(cluster)] = [file_path]

    cluster_list = []
    for key in clusters:
        cluster_list.append(clusters[key])

    threads = args.threads

    process = [None] * args.threads
    all_data_processed = False
    queue = [None] * args.threads

    all_threads_done = False
    thread_done = [False] * args.threads
    count_output = 0
    count_call_of_read_input = 0
    computed_pairs = 0

    matrixFileHandlerObjects_list = [None] * len(cluster_list)

    while not all_data_processed or not all_threads_done:

        for i in range(args.threads):
            if queue[i] is None and not all_data_processed:
                if count_call_of_read_input < len(cluster_list):
                    queue[i] = Queue()
                    process[i] = Process(target=compute_consensus_matrix, kwargs=dict(
                                        pMatrixName = matrices_name,
                                        pMatricesList= cluster_list[count_call_of_read_input], 
                                        pAppend = i > 0, 
                                        pCounter = count_call_of_read_input,
                                        pQueue=queue[i]
                        )
                    )
                    process[i].start()
                    thread_done[i] = False
                else:
                    all_data_processed = True
                    thread_done[i] = True
                count_call_of_read_input += 1
            elif queue[i] is not None and not queue[i].empty():
                matrixFileHandlerObjects, counter_ = queue[i].get()
                matrixFileHandlerObjects_list[counter_] = matrixFileHandlerObjects
                queue[i] = None
                process[i].join()
                process[i].terminate()
                process[i] = None
                thread_done[i] = True
            elif all_data_processed and queue[i] is None:
                thread_done[i] = True
            else:
                time.sleep(1)

        if all_data_processed:
            all_threads_done = True
            for thread in thread_done:
                if not thread:
                    all_threads_done = False

    sum_of_all = []
    for i, matrixFileHandler in enumerate(matrixFileHandlerObjects_list):
        sum_of_all.append(matrixFileHandler.matrixFile.matrix.sum())

    argmin = np.argmin(sum_of_all)

    for i, matrixFileHandler in enumerate(matrixFileHandlerObjects_list):
        # for i, hic_matrix in enumerate(hic_matrix_list):
        matrixFileHandler.matrixFile.matrix.data = matrixFileHandler.matrixFile.matrix.data.astype(np.float32)
        if i != argmin:
            mask = np.isnan(matrixFileHandler.matrixFile.matrix.data)
            matrixFileHandler.matrixFile.matrix.data[mask] = 0

            mask = np.isinf(matrixFileHandler.matrixFile.matrix.data)
            matrixFileHandler.matrixFile.matrix.data[mask] = 0
            adjust_factor = sum_of_all[i] / sum_of_all[argmin]
            matrixFileHandler.matrixFile.matrix.data /= adjust_factor
            mask = np.isnan(matrixFileHandler.matrixFile.matrix.data)

        mask = np.isnan(matrixFileHandler.matrixFile.matrix.data)
        matrixFileHandler.matrixFile.matrix.data[mask] = 0

        mask = np.isinf(matrixFileHandler.matrixFile.matrix.data)
        matrixFileHandler.matrixFile.matrix.data[mask] = 0
        matrixFileHandler.matrixFile.matrix.eliminate_zeros()
        matrixFileHandler.save(args.outFileName + '::/' + 'consensus_matrix_cluster_' + str(i), pSymmetric=True, pApplyCorrection=False)
    # for i, matrixFileHandler in enumerate(matrixFileHandlerObjects_list):
        
