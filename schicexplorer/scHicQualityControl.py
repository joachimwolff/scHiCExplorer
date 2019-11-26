import argparse
from multiprocessing import Process, Queue
import time
import os
import logging
log = logging.getLogger(__name__)

import cooler
from hicmatrix import HiCMatrix as hm

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from schicexplorer._version import __version__


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description=''
    )

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--matrix', '-m',
                                help='The single cell Hi-C interaction matrices to investigate for QC. Needs to be in mcool format',
                                metavar='mcool scHi-C matrix',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--outputMcool', '-o',
                           help='Mcool matrix which contains only the filtered matrices',
                           default='filtered_matrices.mcool')
    parserOpt.add_argument('--minimumReadCoverage',
                           help='Remove all samples with a lower read coverage as this value.',
                           required=False,
                           default=1000000,
                           type=int)
    parserOpt.add_argument('--minimumDensity',
                           help='Remove all samples with a lower read coverage as this value.',
                           required=False,
                           default=0.001,
                           type=float)
    parserOpt.add_argument('--maximumRegionToConsider',
                           help='To compute the density, consider only this genomic distance around the diagonal.',
                           required=False,
                           default=30000000,
                           type=int)
    parserOpt.add_argument('--chromosomes', '-c',
                           nargs='+',
                           help='List of chromosomes that a cell needs to have to be not deleted. However, other chromosomes/contigs and scaffolds which may exist are not deleted. Use scHicAdjustMatrix for this.')
    parserOpt.add_argument('--outFileNameDensity', '-od',
                           help='File name of the density histogram',
                           required=False,
                           default='density.png')
    parserOpt.add_argument('--outFileNameReadCoverage', '-or',
                           help='File name of the read coverage',
                           required=False,
                           default='readCoverage.png')
    parserOpt.add_argument('--outFileNameQCReport', '-oqc',
                           help='File name of the quality report',
                           required=False,
                           default='qc_report.txt')
    parserOpt.add_argument('--dpi', '-d',
                           help='The dpi of the plot.',
                           required=False,
                           default=300,
                           type=int)
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module.',
                           required=False,
                           default=4,
                           type=int)
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def compute_read_coverage_sparsity(pMatrixName, pMatricesList, pXDimension, pMaximumRegionToConsider, pQueue):
    read_coverage = []
    sparsity = []
    for i, matrix in enumerate(pMatricesList):
        hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName + '::' + matrix)

        _matrix = hic_ma.matrix
        max_distance = pMaximumRegionToConsider // hic_ma.getBinSize()

        read_coverage.append(_matrix.data.sum())

        instances, features = _matrix.nonzero()
        distances = np.absolute(instances - features)
        mask = distances >= max_distance
        sparsity_length = len(_matrix.data[mask])

        sparsity.append(sparsity_length / (_matrix.shape[0] * max_distance))

    pQueue.put([read_coverage, sparsity])


def compute_contains_all_chromosomes(pMatrixName, pMatricesList, pChromosomes, pQueue):

    keep_matrices_chromosome_names = []
    for i, matrix in enumerate(pMatricesList):
        ma = hm.hiCMatrix(pMatrixName + '::' + matrix)
        if pChromosomes is None:
            pChromosomes = list(ma.chrBinBoundaries)
        try:
            ma.keepOnlyTheseChr(pChromosomes)
            keep_matrices_chromosome_names.append(1)
        except Exception:
            keep_matrices_chromosome_names.append(0)
    pQueue.put(keep_matrices_chromosome_names)


def main(args=None):

    args = parse_arguments().parse_args(args)

    matrices_name = args.matrix
    threads = args.threads
    matrices_list = cooler.fileops.list_coolers(matrices_name)
    all_samples_number = len(matrices_list)

    #####################################################
    # Detect broken chromosomes and remove these matrices
    #####################################################
    keep_matrices_thread = [None] * threads
    all_data_collected = False
    thread_done = [False] * threads
    length_index = [None] * threads
    length_index[0] = 0
    matricesPerThread = len(matrices_list) // threads
    queue = [None] * threads
    process = [None] * threads
    for i in range(threads):

        if i < threads - 1:
            matrices_name_list = matrices_list[i * matricesPerThread:(i + 1) * matricesPerThread]
            length_index[i + 1] = length_index[i] + len(matrices_name_list)
        else:
            matrices_name_list = matrices_list[i * matricesPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=compute_contains_all_chromosomes, kwargs=dict(
            pMatrixName=matrices_name,
            pMatricesList=matrices_name_list,
            pChromosomes=args.chromosomes,
            pQueue=queue[i]
        )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                worker_result = queue[i].get()
                keep_matrices_thread[i] = worker_result
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

    keep_matrices_chromosome_names = np.array([item for sublist in keep_matrices_thread for item in sublist], dtype=bool)

    matrices_name_chromosome_names = np.array(matrices_list)
    matrices_list = matrices_name_chromosome_names[keep_matrices_chromosome_names]

    matrices_remove = matrices_name_chromosome_names[~keep_matrices_chromosome_names]

    #######################################

    read_coverage = [None] * threads
    sparsity = [None] * threads

    all_data_collected = False
    thread_done = [False] * threads
    length_index = [None] * threads
    length_index[0] = 0
    matricesPerThread = len(matrices_list) // threads
    queue = [None] * threads
    process = [None] * threads
    for i in range(threads):

        if i < threads - 1:
            matrices_name_list = matrices_list[i * matricesPerThread:(i + 1) * matricesPerThread]
            length_index[i + 1] = length_index[i] + len(matrices_name_list)
        else:
            matrices_name_list = matrices_list[i * matricesPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=compute_read_coverage_sparsity, kwargs=dict(
            pMatrixName=matrices_name,
            pMatricesList=matrices_name_list,
            pXDimension=len(matrices_list),
            pMaximumRegionToConsider=args.maximumRegionToConsider,
            pQueue=queue[i]
        )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                worker_result = queue[i].get()
                read_coverage[i] = worker_result[0]
                sparsity[i] = worker_result[1]

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

    read_coverage = np.array([item for sublist in read_coverage for item in sublist])
    sparsity = np.array([item for sublist in sparsity for item in sublist])

    plt.hist(read_coverage, bins=100)
    plt.suptitle('Read coverage of {}'.format(os.path.basename(args.matrix)), fontsize=12)
    plt.title('Matrices with a read coverage < {} are removed.'.format(args.minimumReadCoverage), fontsize=10)
    plt.grid(True)
    plt.axvline(args.minimumReadCoverage, color='r', linestyle='dashed', linewidth=1)
    plt.xlabel('Read coverage')
    plt.ylabel('Frequency')
    plt.savefig(args.outFileNameReadCoverage, dpi=args.dpi)
    plt.close()

    plt.hist(sparsity, bins=100)
    plt.suptitle('Density of {}'.format(os.path.basename(args.matrix)), fontsize=12)
    plt.title('Matrices with a read coverage < {} are removed.'.format(args.minimumReadCoverage), fontsize=10)
    plt.grid(True)
    plt.xlabel('Density')
    plt.ylabel('Frequency')

    plt.axvline(args.minimumDensity, color='r', linestyle='dashed', linewidth=1)

    plt.savefig(args.outFileNameDensity, dpi=args.dpi)
    plt.close()

    mask_read_coverage = read_coverage >= args.minimumReadCoverage
    sum_read_coverage = np.sum(~mask_read_coverage)
    mask_sparsity = sparsity >= args.minimumDensity
    sum_sparsity = np.sum(~mask_sparsity)

    mask = np.logical_or(mask_read_coverage, mask_sparsity)

    matrices_list_filtered = np.array(matrices_list)[mask]

    np.savetxt('accepted_matrices.txt', matrices_list_filtered, fmt="%s")
    np.savetxt('rejected_matrices.txt', np.array(matrices_list)[~mask], fmt="%s")

    if os.path.exists(args.outputMcool):
        os.remove(args.outputMcool)
    for matrix in matrices_list_filtered:

        cooler.fileops.cp(args.matrix + '::' + matrix, args.outputMcool + '::' + matrix)

    ##################
    # Create QC report
    ##################

    header = '# QC report for single-cell Hi-C data generated by scHiCExplorer ' + __version__ + '\n'

    matrix_statistics = 'scHi-C sample contained {} cells:\n'.format(all_samples_number)
    matrices_bad_chromosomes = 'Number of removed matrices containing bad chromosomes {}\n'.format(len(matrices_remove))

    matrices_low_read_coverage = 'Number of removed matrices due to low read coverage (< {}): {}\n'.format(args.minimumReadCoverage, sum_read_coverage)
    matrices_too_sparse = 'Number of removed matrices due to too many zero bins (< {} density, within {} relative genomic distance): {}\n'.format(args.minimumDensity, args.maximumRegionToConsider, sum_sparsity)

    with open(args.outFileNameQCReport, 'w') as file:
        file.write(header)
        file.write(matrix_statistics)
        file.write(matrices_bad_chromosomes)
        file.write(matrices_low_read_coverage)
        file.write(matrices_too_sparse)
