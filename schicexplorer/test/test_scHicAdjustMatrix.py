import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp

from schicexplorer import scHicAdjustMatrix
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test-data/")

import cooler
import numpy.testing as nt
from hicmatrix import HiCMatrix as hm

def test_adjust_matrices_keep():

    outfile = NamedTemporaryFile(suffix='.mcool', delete=False)
    outfile.close()

    chromosomes_to_keep = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX"
    args = "--matrix {} --outFileName {} --action {} --chromosomes {} -t {}".format(ROOT + 'test_matrix.mcool',
                                outfile.name, 'keep', chromosomes_to_keep, 1).split()
    scHicAdjustMatrix.main(args)

    test_data_matrix = ROOT + 'scHicAdjustMatrix/test_matrix_adjusted.mcool'
    matrices_list_test_data = cooler.fileops.list_coolers(test_data_matrix)
    matrices_list_created = cooler.fileops.list_coolers(outfile.name)

    matrices_list_test_data = sorted(matrices_list_test_data)
    matrices_list_created = sorted(matrices_list_created)

    chromosomes_to_keep = sorted(chromosomes_to_keep.split(' '))
    for test_matrix, created_matrix in zip(matrices_list_test_data, matrices_list_created):
        test = hm.hiCMatrix(test_data_matrix + '::' + test_matrix)
        created = hm.hiCMatrix(outfile.name + '::' + created_matrix)
        nt.assert_almost_equal(test.matrix.data, created.matrix.data, decimal=5)
        nt.assert_equal(test.cut_intervals, created.cut_intervals)
        chromosomes_list_test = sorted(cooler.Cooler(test_data_matrix + '::' + test_matrix).chromnames)
        chromosomes_list_created = sorted(cooler.Cooler(outfile.name + '::' + created_matrix).chromnames)
        assert chromosomes_list_test == chromosomes_list_created
        assert chromosomes_list_created == chromosomes_to_keep

        chromosomes_list_test_original = sorted(cooler.Cooler(ROOT + 'test_matrix.mcool' + '::' + test_matrix).chromnames)
        assert chromosomes_list_created != chromosomes_list_test_original
    os.unlink(outfile.name)

def test_adjust_matrices_remove():
    
    outfile = NamedTemporaryFile(suffix='.mcool', delete=False)
    outfile.close()

    chromosomes_to_remove = "chr1 chr2"
    args = "--matrix {} --outFileName {} --action {} --chromosomes {} -t {}".format(ROOT + 'test_matrix.mcool',
                                outfile.name, 'remove', chromosomes_to_remove, 2).split()
    scHicAdjustMatrix.main(args)

    test_data_matrix = ROOT + 'scHicAdjustMatrix/test_matrix_adjusted_remove.mcool'
    matrices_list_test_data = cooler.fileops.list_coolers(test_data_matrix)
    matrices_list_created = cooler.fileops.list_coolers(outfile.name)

    matrices_list_test_data = sorted(matrices_list_test_data)
    matrices_list_created = sorted(matrices_list_created)

    chromosomes_to_remove = sorted(chromosomes_to_remove.split(' '))
    for test_matrix, created_matrix in zip(matrices_list_test_data, matrices_list_created):
        test = hm.hiCMatrix(test_data_matrix + '::' + test_matrix)
        created = hm.hiCMatrix(outfile.name + '::' + created_matrix)
        nt.assert_almost_equal(test.matrix.data, created.matrix.data, decimal=5)
        nt.assert_equal(test.cut_intervals, created.cut_intervals)
        chromosomes_list_test = sorted(cooler.Cooler(test_data_matrix + '::' + test_matrix).chromnames)
        chromosomes_list_created = sorted(cooler.Cooler(outfile.name + '::' + created_matrix).chromnames)
        assert chromosomes_list_test == chromosomes_list_created
        assert chromosomes_to_remove[0] not in chromosomes_list_created
        assert chromosomes_to_remove[1] not in chromosomes_list_created

        chromosomes_list_test_original = sorted(cooler.Cooler(ROOT + 'test_matrix.mcool' + '::' + test_matrix).chromnames)
        assert chromosomes_list_created != chromosomes_list_test_original

    os.unlink(outfile.name)