import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
import numpy.testing as nt

from tempfile import NamedTemporaryFile, mkdtemp

from hicmatrix import HiCMatrix as hm


from schicexplorer import scHicCreateBulkMatrix
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test-data/")


def test_correct_matrices():
    outfile = NamedTemporaryFile(suffix='.scool', delete=False)

    outfile.close()
    args = "--matrix {} --outFileName {} -t {} ".format(ROOT + 'test_matrix.scool',
                                                        outfile.name, 1).split()
    scHicCreateBulkMatrix.main(args)

    test_data_matrix = ROOT + 'scHicCreateBulkMatrix/test_matrix_bulk.cool'

    test = hm.hiCMatrix(test_data_matrix)
    created = hm.hiCMatrix(outfile.name)
    nt.assert_almost_equal(test.matrix.data, created.matrix.data, decimal=5)
    nt.assert_equal(test.cut_intervals, created.cut_intervals)

    os.unlink(outfile.name)


def test_version():
    args = "--version".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicCreateBulkMatrix.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0


def test_help():
    args = "--help".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicCreateBulkMatrix.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0
