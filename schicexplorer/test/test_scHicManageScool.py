import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp
import cooler
import numpy.testing as nt
from hicmatrix import HiCMatrix as hm

from schicexplorer import scHicManageScool
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test-data/scHicManageScool/")


def test_update_matrices():
    outfile = NamedTemporaryFile(suffix='.scool', delete=False)

    outfile.close()
    args = "--matrix {} --outFileName {} -t {} --action update".format(ROOT + 'test_matrix_old.scool',
                                                                       outfile.name, 4, ).split()
    scHicManageScool.main(args)

    test_data_matrix = ROOT + 'test_matrix.scool'
    matrices_list_test_data = cooler.fileops.list_scool_cells(test_data_matrix)
    matrices_list_created = cooler.fileops.list_scool_cells(outfile.name)

    matrices_list_test_data = sorted(matrices_list_test_data)
    matrices_list_created = sorted(matrices_list_created)

    assert len(matrices_list_test_data) == len(matrices_list_created)
    for test_matrix, created_matrix in zip(matrices_list_test_data, matrices_list_created):
        test = hm.hiCMatrix(test_data_matrix + '::' + test_matrix)
        created = hm.hiCMatrix(outfile.name + '::' + created_matrix)
        nt.assert_almost_equal(test.matrix.data, created.matrix.data, decimal=5)
        nt.assert_equal(test.cut_intervals, created.cut_intervals)

    test_data_matrix = ROOT + 'test_matrix_old.scool'
    matrices_list_test_data = cooler.fileops.list_coolers(test_data_matrix)
    matrices_list_created = cooler.fileops.list_scool_cells(outfile.name)

    matrices_list_test_data = sorted(matrices_list_test_data)
    matrices_list_created = sorted(matrices_list_created)

    assert len(matrices_list_test_data) == len(matrices_list_created)
    for test_matrix, created_matrix in zip(matrices_list_test_data, matrices_list_created):
        test = hm.hiCMatrix(test_data_matrix + '::' + test_matrix)
        created = hm.hiCMatrix(outfile.name + '::' + created_matrix)
        nt.assert_almost_equal(test.matrix.data, created.matrix.data, decimal=5)
        nt.assert_equal(test.cut_intervals, created.cut_intervals)

    os.unlink(outfile.name)


def test_extract_matrix():
    output_folder = mkdtemp(prefix="output_")
    args = "--matrix {} --outFileName {} -t {} --action extractToCool -cl {}".format(ROOT + 'test_matrix.scool',
                                                                                     output_folder, 4, ROOT + 'to_extract.txt').split()
    scHicManageScool.main(args)

    matrices_list_test_data = ['Diploid_1_CGTACTAG_CTAAGCCT_R1fastqgz', 'Diploid_1_CGTACTAG_CTCTCTAT_R1fastqgz', 'Diploid_1_CGTACTAG_GTAAGGAG_R1fastqgz']

    for test_matrix in matrices_list_test_data:
        test = hm.hiCMatrix(ROOT + 'test_matrix.scool::' + '/cells/' + test_matrix)
        created = hm.hiCMatrix(output_folder + '/' + test_matrix + '.cool')
        nt.assert_almost_equal(test.matrix.data, created.matrix.data, decimal=5)
        nt.assert_equal(test.cut_intervals, created.cut_intervals)


def test_version():
    args = "--version".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicManageScool.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0


def test_help():
    args = "--help".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicManageScool.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0
