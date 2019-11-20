import argparse
from multiprocessing import Process, Queue
import time

import logging
log = logging.getLogger(__name__)

import cooler
from hicmatrix import HiCMatrix as hm

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


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

    parserRequired.add_argument('--outputMcool', '-o',
                                help='Mcool matrix which contains only the filtered matrices',
                                default='filtered_matrices.mcool')
    parserRequired.add_argument('--minimumReadCoverage',
                                help='Remove all samples with a lower read coverage as this value.',
                                required=False,
                                default=1000000,
                                type=int)
    parserRequired.add_argument('--minimumDensity',
                                help='Remove all samples with a lower read coverage as this value.',
                                required=False,
                                default=0.001,
                                type=float)
    parserRequired.add_argument('--maximumRegionToConsider',
                                help='To compute the density, consider only this genomic distance around the diagonal.',
                                required=False,
                                default=30000000,
                                type=int)
    parserRequired.add_argument('--chromosomes', '-c',
                                nargs='+',
                                help='List of chromosomes that a cell needs to have to be not deleted. However, other chromosomes/contigs and scaffolds which may exist are not deleted. Use scHicAdjustMatrix for this.')
    parserRequired.add_argument('--outFileNameSparsity', '-os',
                                help='File name of the sparsity histogram',
                                required=False,
                                default='sparsity.png')
    parserRequired.add_argument('--outFileNameReadCoverage', '-or',
                                help='File name of the read coverage',
                                required=False,
                                default='readCoverage.png')
    parserRequired.add_argument('--threads', '-t',
                                help='Number of threads. Using the python multiprocessing module.',
                                required=False,
                                default=4,
                                type=int)

    return parser


def compute_read_coverage_sparsity(pMatrixName, pMatricesList, pIndex, pXDimension, pMaximumRegionToConsider, pQueue):
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
    for matrix in pMatricesList:
        ma = hm.hiCMatrix(pMatrixName + '::' + matrix)
        ma.maskBins(ma.nan_bins)
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
    if args.chromosomes:
        keep_matrices_thread = [None] * threads

        all_data_collected = False
        thread_done = [False] * threads
        log.debug('matrix read, starting processing')
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
        np.savetxt('removed_matrices_chromsomes.txt', matrices_remove, fmt="%s")
        log.debug('matrices_remove {} '.format(len(matrices_remove)))
        log.debug('all matrices {} '.format(len(keep_matrices_chromosome_names)))
        log.debug('matrices_list {} '.format(len(matrices_list)))

    #######################################

    read_coverage = [None] * threads
    sparsity = [None] * threads

    all_data_collected = False
    thread_done = [False] * threads
    log.debug('matrix read, starting processing')
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
            pIndex=length_index[i],
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
    plt.savefig(args.outFileNameReadCoverage, dpi=300)
    plt.close()
    plt.hist(sparsity, bins=100)
    plt.savefig(args.outFileNameSparsity, dpi=300)
    plt.close()

    mask_read_coverage = read_coverage >= args.minimumReadCoverage

    mask_sparsity = sparsity >= args.minimumDensity
    mask = np.logical_or(mask_read_coverage, mask_sparsity)
    matrices_list_filtered = np.array(matrices_list)[mask]

    log.debug('len(matrices_list_filtered) {}'.format(len(matrices_list_filtered)))
    log.debug('matrices_list_filtered {}'.format(matrices_list_filtered))
    np.savetxt('accepted_matrices.txt', matrices_list_filtered, fmt="%s")
    np.savetxt('rejected_matrices.txt', np.array(matrices_list)[~mask], fmt="%s")

    for matrix in matrices_list_filtered:

        cooler.fileops.cp(args.matrix + '::' + matrix, args.outputMcool + '::' + matrix)

    ##################
    # Create QC report
    ##################

    # header = '# QC report for single-cell Hi-C data generated by scHiCExplorer ' + __version__ + '\n'

    # matrix_statistics = 'cHi-C sample contained {} cells:'.format(all_samples_number)
    # matrices_bad_chromosomes = 'Number of removed matrices containing bad chromosomes {}\n'.format(len(matrices_removed))

    # matrices_low_read_coverage = 'Number of removed matrices due to low read coverage (< {}): {}\n'.format(args.minimumReadCoverage, len())
    # matrices_too_sparse = 'Number of removed matrices due to too many zero bins (< {} density): {}\n'.format(args.minimumDensity, len())
