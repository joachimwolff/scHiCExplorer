import gzip
import cooler
import re
import os
import time
from multiprocessing import Process, Queue
import numpy as np
import pandas as pd
import logging
log = logging.getLogger(__name__)
from hicmatrix import HiCMatrix as hm
from schicexplorer._version import __version__

from scipy.sparse import csr_matrix
def opener(filename):
    """
    Determines if a file is compressed or not
    """

    f = open(filename, 'rb')
    # print("gzip or not?", f.read(2))

    if f.read(2) == b'\x1f\x8b':
        f.seek(0)
        return gzip.GzipFile(fileobj=f)
    else:
        f.seek(0)
        return f

def cell_name_list(pScoolUri):

    if cooler.fileops.is_scool_file(pScoolUri):
        matrices_list = cooler.fileops.list_scool_cells(pScoolUri)
        r = re.compile('/cell/*')
        # if '/' in matrices_list and any(r.match(line) for line in matrices_list):
        #     matrices_list.remove('/')
        return matrices_list
    else:
        try:
            matrices_list = cooler.fileops.list_coolers(pScoolUri)

            # old and non-standard scool format stored all cells in root
            # no '/' in matrices_list and no '/cell/*'
            log.warning('Please update the scool file to the new file format standard!')
            if not '/' in matrices_list:
                return matrices_list
            # new standard scool format, all cells are stored under '/cell/'
           
            raise Exception('Wrong data format. Please use a scool file.')
            # ['/', '/cells/cell1', '/cells/cell2', '/cells/cell3']

        except Exception:
            raise Exception('Wrong data format. Please use a scool file.')
            exit(1)

def open_and_store_matrix(pMatrixName, pMatricesList, pIndex, pXDimension, pChromosomes, pIntraChromosomalContactsOnly, pChromosomeIndices, pQueue):
    neighborhood_matrix = None
    time_load = 0.0
    time_all = 0.0
    time_csr_create = 0.0
    time_add = 0.0
    features = []
    data = []
    features_length = []
    max_shape = 0
    index_datatype = np.int64
    valid_matrix_list = []
    for i, matrix in enumerate(pMatricesList):
        time_start_load = time.time()
        time_start_all = time.time() 

        cooler_obj = cooler.Cooler(pMatrixName + '::' + matrix)
        shape = cooler_obj.shape
        # pixels = cooler_obj.pixels()[:]
        
        chromosome_indices = None
        mask = None
        chromosome_dataframes_list = []

        # pChromosomes 
        if pChromosomes is None:
            pChromosomes = cooler_obj.chromnames
        for chromosome in pChromosomes:
            
            pixels_chromosome = cooler_obj.pixels().fetch(chromosome)
            # get pixels from one chromosome, but only intra chromosomal contacts       
            # mask_chromosome = pixels['bin1_id'].apply(lambda x: x in pChromosomeIndices[chromosome]) & pixels['bin2_id'].apply(lambda x: x in pChromosomeIndices[chromosome])
            
            # per defintion the bin1_id should belong only to the fetched chromosome, therefore only bin2_id needs to be cleaned
            if pIntraChromosomalContactsOnly:
                mask_chromosome = pixels_chromosome['bin2_id'].apply(lambda x: x in pChromosomeIndices[chromosome])
                chromosome_dataframes_list.append(pixels_chromosome[mask_chromosome])
            else:
                chromosome_dataframes_list.append(pixels_chromosome)

            # if mask is None:
            #     mask = mask_chromosome
            # else:
            #     mask = mask | mask_chromosome
        
        pixels_chromosome = pd.concat(chromosome_dataframes_list)

        if max_shape < shape[0]:
            max_shape = shape[0]
        
        _matrix = [None, None, None]
        if 'bin1_id' in pixels_chromosome.columns and 'bin2_id' in pixels_chromosome.columns  and 'count' in pixels_chromosome.columns :
            _matrix[0] = pixels_chromosome['bin1_id'].values
            _matrix[1] = pixels_chromosome['bin2_id'].values
            _matrix[2] = pixels_chromosome['count'].values
            if len(_matrix[2]) == 0:
                valid_matrix_list.append(False)
                continue
            valid_matrix_list.append(True)
            _matrix[0] = _matrix[0].astype(index_datatype)
            _matrix[1] = _matrix[1].astype(index_datatype)

            _matrix[0] *= np.int64(shape[0]) # matrix[0] are the instance ids, matrix[3] is the shape
            _matrix[0] += _matrix[1] # matrix[3] is the shape, matrix[1] are the feature ids
            features.extend(_matrix[0])
            _matrix[1] = None

            data.extend(_matrix[2])
            features_length.append(len(_matrix[2]))


                # pixels_chromosome['bin1_id'] = pixels_chromosome['bin1_id'] - chromosome_indices[0]
                # pixels_chromosome['bin2_id'] = pixels_chromosome['bin2_id'] - chromosome_indices[0]
            
        # else:
        #     if pChromosomes is not None and len(pChromosomes) == 1:
        #         hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName + '::' + matrix, pChrnameList=pChromosomes, pNoIntervalTree=True, pUpperTriangleOnly=True, pLoadMatrixOnly=True, pRestoreMaskedBins=False)
        #     else:
        #         if not pChromosomes:
        #             hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName + '::' + matrix, pNoIntervalTree=True, pUpperTriangleOnly=True, pLoadMatrixOnly=True, pRestoreMaskedBins=False)
        #         else:
        #             hic_ma = hm.hiCMatrix(pMatrixFile=pMatrixName + '::' + matrix, pNoIntervalTree=False, pUpperTriangleOnly=True, pLoadMatrixOnly=True, pRestoreMaskedBins=False)
        #         if pChromosomes:
        #             hic_ma.keepOnlyTheseChr(pChromosomes)
        #     _matrix = hic_ma.matrix

        #     if len(_matrix[2]) == 0:
        #         valid_matrix_list.append(False)
        #         continue
        #     valid_matrix_list.append(True)
        #     time_load += time.time() - time_start_load
        #     time_csr_create_start = time.time()

        #     time_csr_create += time.time() - time_csr_create_start 
        #     time_add_start = time.time()
        #     if max_shape < _matrix[3]:
        #         max_shape = _matrix[3]
            

        #     _matrix[0] = _matrix[0].astype(index_datatype)
        #     _matrix[1] = _matrix[1].astype(index_datatype)

        #     _matrix[0] *= np.int64(_matrix[3]) # matrix[0] are the instance ids, matrix[3] is the shape
        #     _matrix[0] += _matrix[1] # matrix[3] is the shape, matrix[1] are the feature ids
        #     features.extend(_matrix[0])
        #     _matrix[1] = None

        #     data.extend(_matrix[2])
        #     features_length.append(len(_matrix[2]))
        #     time_add += time.time() - time_add_start
        #     del _matrix
        #     time_all += time.time() - time_start_all
    

    time_start_tocsr = time.time()
    instances = []
    log.debug('')
    for i, instance_length in enumerate(features_length):
        instances.extend([pIndex + i] * instance_length)
    

    neighborhood_matrix = csr_matrix((data, (instances, features)),(pXDimension, max_shape * max_shape), dtype=np.float)
    # log.debug('time_all {}, time_csr {}, time_add {} time_load {} time_tocsr {}'.format(time_all, time_csr_create, time_add, time_load, time.time() - time_start_tocsr))
    pQueue.put([neighborhood_matrix, valid_matrix_list])

def create_csr_matrix_all_cells(pMatrixName, pThreads, pChromosomes, pOutputFolder, pRawFileName, pIntraChromosomalContactsOnly=None):

    matrices_name = pMatrixName
    threads = pThreads
    matrices_list = cell_name_list(matrices_name)
    # log.debug('len(matrices_list) {}'.format(matrices_list))
    neighborhood_matrix = None
    neighborhood_matrix_threads = [None] * threads
    valid_matrix_list_threads  = [None] * threads

    all_data_collected = False
    thread_done = [False] * threads
    # log.debug('matrix read, starting processing')
    length_index = [None] * threads
    length_index[0] = 0
    matricesPerThread = len(matrices_list) // threads
    queue = [None] * threads
    process = [None] * threads

    chromosome_indices = None
    if pIntraChromosomalContactsOnly:
        cooler_obj = cooler.Cooler(pMatrixName + '::' + matrices_list[0])
        binsDataFrame = cooler_obj.bins()[:]
        chromosome_indices = {}
        for chromosome in cooler_obj.chromnames:
            chromosome_indices[chromosome] = np.array(binsDataFrame.index[binsDataFrame['chrom'] == chromosome].tolist())    


    for i in range(threads):

        if i < threads - 1:
            matrices_name_list = matrices_list[i * matricesPerThread:(i + 1) * matricesPerThread]
            length_index[i + 1] = length_index[i] + len(matrices_name_list)
        else:
            matrices_name_list = matrices_list[i * matricesPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=open_and_store_matrix, kwargs=dict(
            pMatrixName=matrices_name,
            pMatricesList=matrices_name_list,
            pIndex=length_index[i],
            pXDimension=len(matrices_list),
            pChromosomes=pChromosomes,
            pIntraChromosomalContactsOnly = pIntraChromosomalContactsOnly,
            pChromosomeIndices=chromosome_indices,
            pQueue=queue[i]
        )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(threads):
            if queue[i] is not None and not queue[i].empty():
                csr_matrix_worker = queue[i].get()
                neighborhood_matrix_threads, valid_matrix_list_threads[i] = csr_matrix_worker
                if neighborhood_matrix is None:
                    neighborhood_matrix = csr_matrix((len(matrices_list), neighborhood_matrix_threads.shape[1]))
                    neighborhood_matrix += neighborhood_matrix_threads
                else:
                    neighborhood_matrix += neighborhood_matrix_threads
                del neighborhood_matrix_threads
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
    # neighborhood_matrix = neighborhood_matrix[neighborhood_matrix.getnnz(1)>0]
    # valid_matrix_list = [item for sublist in valid_matrix_list_threads for item in sublist]
    
    # valid_matrix_list = np.array(valid_matrix_list)
    # matrices_list = np.array(matrices_list)
    # matrices_list_valid = matrices_list[valid_matrix_list]

    # if len(matrices_list_valid) != len(matrices_list):
    #     with open(pOutputFolder + pRawFileName + '_removed_matrices.txt', 'w') as file:
    #         file.write('# Created by scHiCExplorer.utilities version {}\n'.format(__version__))
    #         file.write('# The following cells have not been considered for clustering because they contained no data.\n\n')
    #         file.write("{}".format('\n'.join(matrices_list[~valid_matrix_list])))
    return neighborhood_matrix, matrices_list