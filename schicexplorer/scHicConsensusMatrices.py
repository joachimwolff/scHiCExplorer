import argparse
from multiprocessing import Process, Queue
import time
import logging
log = logging.getLogger(__name__)

import cooler
import numpy as np

from hicmatrix.lib import MatrixFileHandler
from schicexplorer._version import __version__


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description='scHicConsensusMatrices creates based on the clustered samples one consensus matrix for each cluster. '
        'The consensus matrices are normalized to an equal read coverage level and are stored all in one scool matrix.'
    )

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='The single cell Hi-C interaction matrices to investigate for QC. Needs to be in scool format',
                                metavar='scool scHi-C matrix',
                                required=True)
    parserRequired.add_argument('--clusters', '-c',
                                help='Text file which contains per matrix the associated cluster.',
                                metavar='cluster file',
                                required=True)
    parserRequired.add_argument('--outFileName', '-o',
                                help='File name of the consensus scool matrix.',
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


def compute_consensus_matrix(pMatrixName, pClusterMatricesList, pAppend, pQueue):
    cluster_consensus_matrices_list = []
    for i, cluster in enumerate(pClusterMatricesList):
        consensus_matrix = None
        if i == 0 and pAppend:
            append = False
        else:
            append = True
        for matrix in cluster:
            matrixFileHandlerInput = MatrixFileHandler(pFileType='cool', pMatrixFile=pMatrixName + '::' + matrix)
            _matrix, cut_intervals, nan_bins, \
                distance_counts, correction_factors = matrixFileHandlerInput.load()

            if consensus_matrix is None:
                consensus_matrix = _matrix
            else:
                consensus_matrix += _matrix

        hic2CoolVersion = matrixFileHandlerInput.matrixFile.hic2cool_version
        matrixFileHandlerOutput = MatrixFileHandler(pFileType='cool', pAppend=append, pEnforceInteger=False, pFileWasH5=False, pHic2CoolVersion=hic2CoolVersion)

        matrixFileHandlerOutput.set_matrix_variables(consensus_matrix, cut_intervals, nan_bins,
                                                     correction_factors, distance_counts)
        cluster_consensus_matrices_list.append(matrixFileHandlerOutput)

    pQueue.put(cluster_consensus_matrices_list)


def main(args=None):

    args = parse_arguments().parse_args(args)

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
    if len(cluster_list) < threads:
        threads = len(cluster_list)

    consensus_matrices_threads = [None] * threads
    all_data_collected = False
    thread_done = [False] * threads
    length_index = [None] * threads
    length_index[0] = 0
    clusterPerThread = len(cluster_list) // threads
    queue = [None] * threads
    process = [None] * threads
    for i in range(threads):

        if i < threads - 1:
            cluster_name_list = cluster_list[i * clusterPerThread:(i + 1) * clusterPerThread]
        else:
            cluster_name_list = cluster_list[i * clusterPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=compute_consensus_matrix, kwargs=dict(
            pMatrixName=args.matrix,
            pClusterMatricesList=cluster_name_list,
            pAppend=i == 0,
            pQueue=queue[i]
        )
        )
        process[i].start()

    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                consensus_matrices_threads[i] = queue[i].get()

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

    matrixFileHandlerObjects_list = [item for sublist in consensus_matrices_threads for item in sublist]

    sum_of_all = []
    for i, matrixFileHandler in enumerate(matrixFileHandlerObjects_list):
        sum_of_all.append(matrixFileHandler.matrixFile.matrix.sum())

    argmin = np.argmin(sum_of_all)

    for i, matrixFileHandler in enumerate(matrixFileHandlerObjects_list):
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
