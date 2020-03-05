import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp

from schicexplorer import scHicMergeToMCool
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test-data/")


import cooler
import numpy.testing as nt
from hicmatrix import HiCMatrix as hm


def test_correct_matrices():

    outfile = NamedTemporaryFile(suffix='.mcool', delete=False)
    outfile.close()

    input_matrices = "Diploid_1_CGTACTAG_AAGGAGTA_R1fastqgz.cool Diploid_1_CGTACTAG_ACTGCATA_R1fastqgz.cool Diploid_1_CGTACTAG_CGTCTAAT_R1fastqgz.cool Diploid_1_CGTACTAG_CTAAGCCT_R1fastqgz.cool Diploid_1_CGTACTAG_CTCTCTAT_R1fastqgz.cool Diploid_1_CGTACTAG_GTAAGGAG_R1fastqgz.cool Diploid_1_CGTACTAG_TATCCTCT_R1fastqgz.cool Diploid_1_CGTACTAG_TCTCTCCG_R1fastqgz.cool Diploid_1_TAAGGCGA_AAGGAGTA_R1fastqgz.cool Diploid_1_TAAGGCGA_CGTCTAAT_R1fastqgz.cool Diploid_1_TAAGGCGA_CTAAGCCT_R1fastqgz.cool Diploid_2_AAGAGGCA_AAGGAGTA_R1fastqgz.cool Diploid_2_AAGAGGCA_ACTGCATA_R1fastqgz.cool Diploid_2_AAGAGGCA_CGTCTAAT_R1fastqgz.cool Diploid_2_AAGAGGCA_CTAAGCCT_R1fastqgz.cool Diploid_2_AAGAGGCA_CTCTCTAT_R1fastqgz.cool Diploid_2_AAGAGGCA_GTAAGGAG_R1fastqgz.cool Diploid_2_AAGAGGCA_TATCCTCT_R1fastqgz.cool Diploid_2_AAGAGGCA_TCTCTCCG_R1fastqgz.cool Diploid_2_AGGCAGAA_AAGGAGTA_R1fastqgz.cool"
    input_matrices = input_matrices.split(' ')
    input_matrices_str = ''
    for i in range(len(input_matrices)):
        input_matrices_str += ROOT + 'scHicMergeToMCool/' + input_matrices[i] + ' '

    args = "--matrices {} --outFileName {} ".format(input_matrices_str,
                                                    outfile.name).split()
    scHicMergeToMCool.main(args)

    test_data_matrix = ROOT + 'test_matrix.mcool'
    matrices_list_test_data = cooler.fileops.list_coolers(test_data_matrix)
    matrices_list_created = cooler.fileops.list_coolers(outfile.name)

    matrices_list_test_data = sorted(matrices_list_test_data)
    matrices_list_created = sorted(matrices_list_created)

    for test_matrix, created_matrix in zip(matrices_list_test_data, matrices_list_created):
        test = hm.hiCMatrix(test_data_matrix + '::' + test_matrix)
        created = hm.hiCMatrix(outfile.name + '::' + created_matrix)
        nt.assert_almost_equal(test.matrix.data, created.matrix.data, decimal=5)
        nt.assert_equal(test.cut_intervals, created.cut_intervals)

    os.unlink(outfile.name)


def test_version():
    args = "--version".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicMergeToMCool.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0


def test_help():
    args = "--help".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicMergeToMCool.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0
