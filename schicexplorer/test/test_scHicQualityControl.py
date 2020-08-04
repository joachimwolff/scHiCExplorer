import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp

from schicexplorer import scHicQualityControl
import cooler
import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test-data/")
import numpy.testing as nt
from hicmatrix import HiCMatrix as hm
tolerance = 40


def are_files_equal(file1, file2, delta=2, skip=0):
    equal = True
    if delta:
        mismatches = 0
    with open(file1) as textfile1, open(file2) as textfile2:
        for i, (x, y) in enumerate(zip(textfile1, textfile2)):
            if i < skip:
                continue
            if x != y:
                if delta:
                    mismatches += 1
                    if mismatches > delta:
                        equal = False
                        break
                else:
                    equal = False
                    break
    return equal


def test_plot():
    outfile_density = NamedTemporaryFile(suffix='.png', delete=False)
    outfile_density.close()
    outfile_coverage = NamedTemporaryFile(suffix='.png', delete=False)
    outfile_coverage.close()
    outfile_qc_report = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile_qc_report.close()
    outfile_matrix = NamedTemporaryFile(suffix='.scool', delete=False)
    outfile_matrix.close()
    args = "--matrix {} --outputScool {} -t {} --dpi {} --outFileNameDensity {} \
            --outFileNameReadCoverage {} --outFileNameQCReport {} \
            --minimumReadCoverage {} --minimumDensity {} \
            --maximumRegionToConsider {} --runChromosomeCheck".format(ROOT + 'test_matrix.scool',
                                                                      outfile_matrix.name, 1, 300,
                                                                      outfile_density.name,
                                                                      outfile_coverage.name,
                                                                      outfile_qc_report.name,
                                                                      100000, 0.001, 30000000
                                                                      ).split()
    scHicQualityControl.main(args)

    test_image_density = ROOT + 'scHicQualityControl/density.png'
    res = compare_images(test_image_density, outfile_density.name, tolerance)
    assert res is None, res

    test_image_density = ROOT + 'scHicQualityControl/coverage.png'
    res = compare_images(test_image_density, outfile_coverage.name, tolerance)
    assert res is None, res

    assert are_files_equal(ROOT + "scHicQualityControl/qc_report.txt", outfile_qc_report.name)

    test_data_matrix = ROOT + 'scHicQualityControl/qc_matrix.scool'
    matrices_list_test_data = cooler.fileops.list_scool_cells(test_data_matrix)
    matrices_list_created = cooler.fileops.list_scool_cells(outfile_matrix.name)

    matrices_list_test_data = sorted(matrices_list_test_data)
    matrices_list_created = sorted(matrices_list_created)

    for test_matrix, created_matrix in zip(matrices_list_test_data, matrices_list_created):
        test = hm.hiCMatrix(test_data_matrix + '::' + test_matrix)
        created = hm.hiCMatrix(outfile_matrix.name + '::' + created_matrix)
        nt.assert_almost_equal(test.matrix.data, created.matrix.data, decimal=5)
        nt.assert_equal(test.cut_intervals, created.cut_intervals)


def test_plot_chromosomes():
    outfile_density = NamedTemporaryFile(suffix='.png', delete=False)
    outfile_density.close()
    outfile_coverage = NamedTemporaryFile(suffix='.png', delete=False)
    outfile_coverage.close()
    outfile_qc_report = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile_qc_report.close()
    outfile_matrix = NamedTemporaryFile(suffix='.scool', delete=False)
    outfile_matrix.close()
    args = "--matrix {} --outputScool {} -t {} --dpi {} --outFileNameDensity {} \
            --outFileNameReadCoverage {} --outFileNameQCReport {} \
            --minimumReadCoverage {} --minimumDensity {} \
            --maximumRegionToConsider {} --chromosomes chr1 chr2".format(ROOT + 'test_matrix.scool',
                                                                         outfile_matrix.name, 1, 300,
                                                                         outfile_density.name,
                                                                         outfile_coverage.name,
                                                                         outfile_qc_report.name,
                                                                         100000, 0.001, 30000000
                                                                         ).split()
    scHicQualityControl.main(args)

    test_image_density = ROOT + 'scHicQualityControl/density_chr1_chr2.png'
    res = compare_images(test_image_density, outfile_density.name, tolerance)
    assert res is None, res

    test_image_density = ROOT + 'scHicQualityControl/coverage_chr1_chr2.png'
    res = compare_images(test_image_density, outfile_coverage.name, tolerance)
    assert res is None, res

    assert are_files_equal(ROOT + "scHicQualityControl/qc_report_chr1_chr2.txt", outfile_qc_report.name)

    test_data_matrix = ROOT + 'scHicQualityControl/qc_matrix_chr1_chr2.scool'
    matrices_list_test_data = cooler.fileops.list_scool_cells(test_data_matrix)
    matrices_list_created = cooler.fileops.list_scool_cells(outfile_matrix.name)

    matrices_list_test_data = sorted(matrices_list_test_data)
    matrices_list_created = sorted(matrices_list_created)

    for test_matrix, created_matrix in zip(matrices_list_test_data, matrices_list_created):
        test = hm.hiCMatrix(test_data_matrix + '::' + test_matrix)
        created = hm.hiCMatrix(outfile_matrix.name + '::' + created_matrix)
        nt.assert_almost_equal(test.matrix.data, created.matrix.data, decimal=5)
        nt.assert_equal(test.cut_intervals, created.cut_intervals)


def test_version():
    args = "--version".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicQualityControl.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0


def test_help():
    args = "--help".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicQualityControl.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0
