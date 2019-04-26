# read all matrices
## get non-zeros and flatten it (x*length) + y
# make number of instacnes * dim**2 csr matrix

import argparse
import os
import gzip
import shutil
from multiprocessing import Process, Queue
import time
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

import logging
log = logging.getLogger(__name__)

import cooler
from sparse_neighbors_search import MinHash
from sparse_neighbors_search import MinHashDBSCAN
from sparse_neighbors_search import MinHashSpectralClustering


from hicmatrix import HiCMatrix
from schicexplorer.utilities import opener
from hicmatrix.lib import MatrixFileHandler

from scipy.sparse import csr_matrix
def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        add_help=False
    )

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--matrix', '-m',
                                help='The single cell Hi-C interaction matrices to cluster. Needs to be in mcool format',
                                metavar='mcool scHi-C matrix',
                                required=True)

    parserRequired.add_argument('--outputFolder', '-o',
                                help='Path of folder to save the demultiplexed files',
                                metavar='FOLDER',
                                required=False,
                                default='clusters')

    parserRequired.add_argument('--threads',
                           help='Number of threads. Using the python multiprocessing module.',
                           required=False,
                           default=4,
                           type=int)

    return parser

def mapNeihgborsToMatrixNames():
    pass

def main(args=None):

    args = parse_arguments().parse_args(args)
    # print('hello world from scHiCExplorer')
    if not os.path.exists(args.outputFolder):
        try:
            os.makedirs(args.outputFolder)
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    matrices = cooler.fileops.list_coolers(args.matrix)
    log.debug('matrices {}'.format(matrices))
    neighborhood_matrix = None
    for i, matrix in enumerate(matrices):
        matrixFileHandlerInput = MatrixFileHandler(pFileType='cool', pMatrixFile=args.matrix+ '::' +matrix)
        _matrix, _, _, \
                _, _ = matrixFileHandlerInput.load()
    
        if neighborhood_matrix is None:
            neighborhood_matrix = csr_matrix((len(matrices), _matrix.shape[0] * _matrix.shape[1]), dtype=int)

        instances, features = _matrix.nonzero()

        instances *= _matrix.shape[1]
        instances += features
        features = None
        neighborhood_matrix[i, instances] = 1 #TODO this is so far wrong

    # log.debug('neighborhood_matrix {}'.format(neighborhood_matrix))
    minHash = MinHash(number_of_hash_functions=20, number_of_cores=4)
    minHash.fit(neighborhood_matrix)
    print(minHash.kneighbors_graph(mode='distance'))
    minHashDbscan = MinHashDBSCAN(number_of_hash_functions=20, number_of_cores=4)
    # minHashClustering = MinHashClustering(minHash)
    print(minHashDbscan.fit_predict(neighborhood_matrix))

    minHashSpectralClustering = MinHashSpectralClustering(number_of_hash_functions=20, number_of_cores=4)
    # minHashClustering = MinHashClustering(minHash)
    print(minHashSpectralClustering.fit_predict(neighborhood_matrix))

    # with