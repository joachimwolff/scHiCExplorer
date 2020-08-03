import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp
import cooler
import numpy.testing as nt
from hicmatrix import HiCMatrix as hm

from schicexplorer import scHicConsensusMatrices
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test-data/")


def test_consensus_matrices():
    outfile_cell_names = NamedTemporaryFile(prefix='.txt', delete=False)
    outfile_chromosome_sizes = NamedTemporaryFile(prefix='.txt', delete=False)
    output_folder = mkdtemp(prefix="output_")

    outfile.close()
    args = "--matrix {} --outFileName {} -t {} -c {}".format(ROOT + 'test_matrix.scool',
                                                             outfile.name, 4, ROOT + 'scHicConsensusMatrices/cluster_kmeans.txt').split()
    scHicConsensusMatrices.main(args)

    test_data_matrix = ROOT + 'scHicConsensusMatrices/consensus_matrix.scool'
    matrices_list_test_data = cooler.fileops.list_coolers(test_data_matrix)
    matrices_list_created = cooler.fileops.list_coolers(outfile.name)

    matrices_list_test_data = sorted(matrices_list_test_data)
    matrices_list_created = sorted(matrices_list_created)

    assert len(matrices_list_test_data) == len(matrices_list_created)
    for test_matrix, created_matrix in zip(matrices_list_test_data, matrices_list_created):
        test = hm.hiCMatrix(test_data_matrix + '::' + test_matrix)
        created = hm.hiCMatrix(outfile.name + '::' + created_matrix)
        nt.assert_almost_equal(test.matrix.data, created.matrix.data, decimal=5)
        nt.assert_equal(test.cut_intervals, created.cut_intervals)

    os.unlink(outfile.name)


def test_version():
    args = "--version".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicConsensusMatrices.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0


def test_help():
    args = "--help".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicConsensusMatrices.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0
