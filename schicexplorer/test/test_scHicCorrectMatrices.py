import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp

import cooler
import numpy.testing as nt
from hicmatrix import HiCMatrix as hm

from schicexplorer import scHicCorrectMatrices
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test-data/")



def test_correct_matrices():
    outfile = NamedTemporaryFile(suffix='.mcool', delete=False)

    outfile.close()
    args = "--matrix {} --outFileName {} -t {} ".format(ROOT + 'test_matrix.mcool',
                                outfile.name, 1).split()
    scHicCorrectMatrices.main(args)

    test_data_matrix = ROOT + 'scHicCorrectMatrices/test_matrix_corrected.mcool'
    matrices_list_test_data = cooler.fileops.list_coolers(test_data_matrix)
    matrices_list_created = cooler.fileops.list_coolers(outfile.name)

    matrices_list_test_data = sorted(matrices_list_test_data)
    matrices_list_created = sorted(matrices_list_created)

    for test_matrix, created_matrix in zip(matrices_list_test_data, matrices_list_created):
        test = hm.hiCMatrix(test_data_matrix + '::' + test_matrix)
        created = hm.hiCMatrix(outfile.name + '::' + created_matrix )
        nt.assert_almost_equal(test.matrix.data, created.matrix.data, decimal=5)
        nt.assert_equal(test.cut_intervals, created.cut_intervals)

    os.unlink(outfile.name)