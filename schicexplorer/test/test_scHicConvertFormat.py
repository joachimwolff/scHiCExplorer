import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp
import cooler
import numpy.testing as nt
from hicmatrix import HiCMatrix as hm

from schicexplorer import scHicConvertFormat
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test-data/")


def are_files_equal(file1, file2, delta=2, skip=0):
    equal = True
    if delta:
        mismatches = 0
    with open(file1) as textfile1, open(file2) as textfile2:
        for i, (x, y) in enumerate(zip(textfile1, textfile2)):
            if i < skip:
                continue
            if x[-10:] == y[-10:]:
                continue
            else:
                if delta:
                    mismatches += 1
                    if mismatches > delta:
                        equal = False
                        break
                else:
                    equal = False
                    break
    return equal


def test_consensus_matrices():
    outfile_cell_names = NamedTemporaryFile(prefix='.txt', delete=False)
    outfile_chromosome_sizes = NamedTemporaryFile(prefix='.txt', delete=False)
    output_folder = mkdtemp(prefix="output_")

    outfile_cell_names.close()
    outfile_chromosome_sizes.close()

    args = "--matrix {} -t {} -of {} -oc {} -os {}".format(ROOT + 'test_matrix.scool', 4,
                                                           output_folder, outfile_cell_names.name, outfile_chromosome_sizes.name).split()
    scHicConvertFormat.main(args)

    assert set(os.listdir(ROOT + "scHicConvertFormat/scHiCluster/cells/")) == set(os.listdir(output_folder + '/cells/'))
    assert are_files_equal(ROOT + 'scHicConvertFormat/scHiCluster/cellNameFile.txt', outfile_cell_names.name, delta=2, skip=0)
    assert are_files_equal(ROOT + 'scHicConvertFormat/scHiCluster/chromosomeSize.txt', outfile_chromosome_sizes.name, delta=2, skip=0)


def test_version():
    args = "--version".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicConvertFormat.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0


def test_help():
    args = "--help".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicConvertFormat.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0
