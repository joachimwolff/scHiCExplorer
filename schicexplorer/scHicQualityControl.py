import argparse
from multiprocessing import Process, Queue
import time
import os
import logging
log = logging.getLogger(__name__)

import cooler
from hicmatrix import HiCMatrix as hm
from hicmatrix.lib import MatrixFileHandler
from datetime import datetime

import numpy as np
from scipy.sparse import csr_matrix
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from schicexplorer._version import __version__
from schicexplorer.utilities import cell_name_list


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description=''
    )

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--matrix', '-m',
                                help='The single cell Hi-C interaction matrices to investigate for QC. Needs to be in scool format',
                                metavar='scool scHi-C matrix',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--outputScool', '-o',
                           help='scool matrix which contains only the filtered matrices',
                           default='filtered_matrices.scool')
    parserOpt.add_argument('--minimumReadCoverage',
                           help='Remove all samples with a lower read coverage as this value.',
                           required=False,
                           default=1000000,
                           type=int)
    parserOpt.add_argument('--minimumDensity',
                           help='Remove all samples with a lower density as this value. The density is given by: number of non-zero interactions / all possible interactions.',
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
    parserOpt.add_argument('--plotOnly',
                           help='Do not create a new matrix, create only the plots.',
                           action='store_true')
    parserOpt.add_argument('--runChromosomeCheck',
                           help='Skip the data integrity check for the chromosomes.',
                           action='store_true')
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
    # matrixFileHandlerList = []

    log.debug('read covarage and sparsity')
    hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName + '::' + pMatricesList[0])
    bin_size = hic_ma.getBinSize()
    shape_x = hic_ma.matrix.shape[0]
    for i, matrix in enumerate(pMatricesList):

        # _matrix = hic_ma.matrix

        matrixFileHandler = MatrixFileHandler(pFileType='cool', pMatrixFile=pMatrixName + '::' + matrix, pLoadMatrixOnly=True)
        _matrix, cut_intervals, nan_bins, \
            distance_counts, correction_factors = matrixFileHandler.load()
        max_distance = pMaximumRegionToConsider // bin_size


        # if (_matrix.data.sum()) < 100000:
        #     log.debug('smaller 100000: {}'.format(_matrix.data.sum()))
        # if np.sum(_matrix.data) < 100000:
        #     log.debug('II smaller 100000: {}'.format(_matrix.data.sum()))
        instances = _matrix[0]
        features = _matrix[1]

        distances = np.absolute(instances - features)
        mask = distances <= max_distance
        sparsity_length = len(_matrix[2][mask])

        sparsity.append(sparsity_length / (shape_x * max_distance))
        
        # only upper half is loaded --> times 2
        read_coverage_sum = _matrix[2].sum()*2 
        # minus the double main diagonal
        mask = distances == 0
        read_coverage_sum -= _matrix[2][mask].sum()
        read_coverage.append(read_coverage_sum)
        # matrixFileHandlerOutput = MatrixFileHandler(pFileType='cool', pMatrixFile=matrix, pAppend=False, pEnforceInteger=False, pFileWasH5=False, pHic2CoolVersion=None)

        # _matrix = csr_matrix((_matrix[2], (_matrix[0], _matrix[1])),(_matrix[3], _matrix[3]), dtype=np.float)
        # matrixFileHandlerOutput.set_matrix_variables(_matrix, cut_intervals, nan_bins,
        #                                              correction_factors, distance_counts)
        # log.debug('\n\nnew read after create foo: {}'.format(np.sum(_matrix.data)))

        # matrixFileHandlerList.append(matrixFileHandlerOutput)

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
    matrices_list = cell_name_list(matrices_name)
    all_samples_number = len(matrices_list)

    if args.runChromosomeCheck:
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
    # log.debug('matrices_remove {}'.format(matrices_remove[:10]))

    #######################################

    read_coverage_thread = [None] * threads
    sparsity_thread = [None] * threads

    all_data_collected = False
    thread_done = [False] * threads
    length_index = [None] * threads
    matrixFileObjects_thread = [None] * threads
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
                read_coverage_thread[i] = worker_result[0]
                sparsity_thread[i] = worker_result[1]
                # matrixFileObjects_thread[i] = worker_result[2]
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

    read_coverage = np.array([item for sublist in read_coverage_thread for item in sublist])
    sparsity = np.array([item for sublist in sparsity_thread for item in sublist])
    # matrixFileObjects = np.array([item for sublist in matrixFileObjects_thread for item in sublist])


    log.debug('read_coverage {}'.format(read_coverage))
    plt.close()
    plt.hist(read_coverage, bins=100)

    plt.suptitle('Read coverage of {}'.format(os.path.basename(args.matrix)), fontsize=12)
    plt.grid(True)
    if args.minimumReadCoverage > 0:
        plt.axvline(args.minimumReadCoverage, color='r', linestyle='dashed', linewidth=1)
        plt.title('Matrices with a read coverage < {} are removed.'.format(args.minimumReadCoverage), fontsize=10)

    plt.xlabel('Read coverage')
    plt.ylabel('Frequency')
    plt.savefig(args.outFileNameReadCoverage, dpi=args.dpi)
    plt.close()

    plt.hist(sparsity, bins=100)
    plt.suptitle('Density of {}'.format(os.path.basename(args.matrix)), fontsize=12)
    if args.minimumDensity > 0:
        plt.title('Matrices with a density < {} are removed.'.format(args.minimumDensity), fontsize=10)
    plt.grid(True)
    plt.xlabel('Density')
    plt.ylabel('Frequency')
    if args.minimumDensity > 0:
        plt.axvline(args.minimumDensity, color='r', linestyle='dashed', linewidth=1)

    plt.savefig(args.outFileNameDensity, dpi=args.dpi)
    plt.close()

    mask_read_coverage = read_coverage >= args.minimumReadCoverage
    mask_sparsity = sparsity >= args.minimumDensity

    log.debug('len(mask_read_coverage) {}'.format(len(mask_read_coverage)))
    log.debug('len(mask_sparsity) {}'.format(len(mask_sparsity)))

    mask = np.logical_and(mask_read_coverage, mask_sparsity)
    # mask = np.logical_or(mask, matrices_remove)
    matrices_list_filtered = np.array(matrices_list)[mask]
    # matrixFile_objects_filtered = np.array(matrixFileObjects)[mask]


    # log.debug('matrixFile_objects_filtered {}'.format(matrixFile_objects_filtered))
    # matrixFileObjects
    sum_read_coverage = np.sum(~mask_read_coverage)
    sum_sparsity = np.sum(~mask_sparsity)

    if not args.plotOnly:
        np.savetxt('accepted_matrices.txt', matrices_list_filtered, fmt="%s")
        np.savetxt('rejected_matrices.txt', np.array(matrices_list)[~mask], fmt="%s")

        if os.path.exists(args.outputScool):
            os.remove(args.outputScool)
        log.debug('matrices_list_filtered {}'.format(matrices_list_filtered[:10]))

        cooler.fileops.cp(args.matrix + '::/bins', args.outputScool + '::/bins')
        cooler.fileops.cp(args.matrix + '::/chroms', args.outputScool + '::/chroms')

        
        # with cooler.util.open_hdf5(path) as h5:
        


        with cooler.util.open_hdf5(args.matrix) as source:
            attributes_dict = {}
            for k, v in source.attrs.items():
                attributes_dict[k] = v

            attributes_dict['ncells'] = len(matrices_list_filtered)
            attributes_dict['creation-date'] = datetime.now().isoformat()
            with h5py.File(args.outputScool, "r+") as f:
                h5 = f['/']
                h5.attrs.update(attributes_dict)

        content_bins_ln =  ['chrom', 'start', 'end']
        for matrix in matrices_list_filtered:
            
            cooler.fileops.cp(args.matrix + '::' + matrix + '/pixels', args.outputScool + '::' + matrix + '/pixels')
            cooler.fileops.cp(args.matrix + '::' + matrix + '/indexes', args.outputScool + '::' + matrix+ '/indexes')
            cooler.fileops.ln(args.outputScool + '::' + '/chroms', args.outputScool + '::' + matrix+ '/chroms')
            cooler.fileops.ln(args.outputScool + '::' + '/bins/chrom', args.outputScool + '::' + matrix+ '/bins/chrom')
            cooler.fileops.ln(args.outputScool + '::' + '/bins/start', args.outputScool + '::' + matrix+ '/bins/start')
            cooler.fileops.ln(args.outputScool + '::' + '/bins/end', args.outputScool + '::' + matrix+ '/bins/end')

            group_dataset_list = cooler.fileops.ls(args.matrix + '::' + matrix + '/bins/')
            for datatype in group_dataset_list:
                last_element = datatype.split('/')[-1]
                if not (last_element) in content_bins_ln and last_element != '':
                    cooler.fileops.cp(args.matrix + '::' + matrix + '/bins/'+last_element, args.outputScool + '::' + matrix+ '/bins/'+last_element)


            with cooler.util.open_hdf5(args.matrix) as source:#, cooler.util.open_hdf5(args.outputScool + '::' + matrix) as destination:
                # attributes_dict = source.attrs

                attributes_dict = {}
                for k, v in source[matrix].attrs.items():
                    attributes_dict[k] = v
                with h5py.File(args.outputScool, "r+") as f:
                    h5 = f[matrix]
                    h5.attrs.update(attributes_dict)
                # attributes_dict['ncells'] = len(matrices_list_filtered)
                # attributes_dict['creation-date'] = datetime.now().isoformat()
                # destination.attrs.update(attributes_dict)


        # matrixFileHandler = MatrixFileHandler(pFileType='scool')
        # matrixFileHandler.matrixFile.coolObjectsList = matrixFile_objects_filtered
        # matrixFileHandler.save(args.outputScool, pSymmetric=True, pApplyCorrection=False)
  


    ##################
    # Create QC report
    ##################

    header = '# QC report for single-cell Hi-C data generated by scHiCExplorer ' + __version__ + '\n'

    matrix_statistics = 'scHi-C sample contained {} cells:\n'.format(all_samples_number)
    if args.runChromosomeCheck:
        matrices_bad_chromosomes = 'Number of removed matrices containing bad chromosomes {}\n'.format(len(matrices_remove))

    matrices_low_read_coverage = 'Number of removed matrices due to low read coverage (< {}): {}\n'.format(args.minimumReadCoverage, sum_read_coverage)
    matrices_too_sparse = 'Number of removed matrices due to too many zero bins (< {} density, within {} relative genomic distance): {}\n'.format(args.minimumDensity, args.maximumRegionToConsider, sum_sparsity)

    matrix_qc = '{} samples passed the quality control. Please consider matrices with a low read coverage may be the matrices with a low density and overlap therefore.'.format(len(matrices_list_filtered))

    with open(args.outFileNameQCReport, 'w') as file:
        file.write(header)
        file.write(matrix_statistics)
        if args.runChromosomeCheck:
            file.write(matrices_bad_chromosomes)
        file.write(matrices_low_read_coverage)
        file.write(matrices_too_sparse)
        file.write(matrix_qc)
