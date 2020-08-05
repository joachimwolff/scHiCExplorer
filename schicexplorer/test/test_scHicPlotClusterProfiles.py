import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp
import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from schicexplorer import scHicPlotClusterProfiles
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test-data/")

tolerance = 60


def test_plot_svl():
    outfile = NamedTemporaryFile(suffix='.png', delete=False)

    outfile.close()
    args = "--matrix {} -c {} --outFileName {} -t {} --maximalDistance {} --distanceShortRange {} --distanceLongRange {} --orderBy {} --dpi {}".format(ROOT + 'test_matrix.scool',
                                                                                                                                                       ROOT + 'scHicPlotClusterProfiles/cluster_kmeans.txt', outfile.name,
                                                                                                                                                       1, 50000000, 2000000, 12000000,
                                                                                                                                                       'svl', 300

                                                                                                                                                       ).split()
    scHicPlotClusterProfiles.main(args)
    test_image_path = ROOT + 'scHicPlotClusterProfiles/plot.png'
    res = compare_images(test_image_path, outfile.name, tolerance)
    assert res is None, res

    # os.unlink(outfile.name)


def test_plot_use_defaults():
    outfile = NamedTemporaryFile(suffix='.png', delete=False)

    outfile.close()
    args = "--matrix {} -c {} --outFileName {} -t {} ".format(ROOT + 'test_matrix.scool',
                                                              ROOT + 'scHicPlotClusterProfiles/cluster_kmeans.txt', outfile.name,
                                                              1
                                                              ).split()
    scHicPlotClusterProfiles.main(args)
    test_image_path = ROOT + 'scHicPlotClusterProfiles/plot.png'
    res = compare_images(test_image_path, outfile.name, tolerance)
    assert res is None, res

    # os.unlink(outfile.name)


def test_plot_chromosomes():
    outfile = NamedTemporaryFile(suffix='.png', delete=False)

    outfile.close()
    args = "--matrix {} -c {} --outFileName {} -t {} --maximalDistance {} --distanceShortRange {} \
            --distanceLongRange {} --orderBy {} --dpi {} --chromosomes {}"\
                .format(ROOT + 'test_matrix.scool',
                        ROOT + 'scHicPlotClusterProfiles/cluster_kmeans.txt', outfile.name,
                        1, 50000000, 2000000, 12000000,
                        'svl', 300, 'chr1 chr2'
                        ).split()
    scHicPlotClusterProfiles.main(args)
    test_image_path = ROOT + 'scHicPlotClusterProfiles/plot_chr1_chr2.png'
    res = compare_images(test_image_path, outfile.name, tolerance)
    assert res is None, res

    # os.unlink(outfile.name)


def test_plot_orderByFile():
    outfile = NamedTemporaryFile(suffix='.png', delete=False)

    outfile.close()
    args = "--matrix {} -c {} --outFileName {} -t {} --maximalDistance {} --distanceShortRange {} \
            --distanceLongRange {} --orderBy {} --dpi {}"\
                .format(ROOT + 'test_matrix.scool',
                        ROOT + 'scHicPlotClusterProfiles/cluster_kmeans.txt', outfile.name,
                        1, 50000000, 2000000, 12000000,
                        'orderByFile', 300
                        ).split()
    scHicPlotClusterProfiles.main(args)
    test_image_path = ROOT + 'scHicPlotClusterProfiles/plot_orderByFile.png'
    res = compare_images(test_image_path, outfile.name, tolerance)
    assert res is None, res

    # os.unlink(outfile.name)


def test_plot_maxDistance():
    outfile = NamedTemporaryFile(suffix='.png', delete=False)

    outfile.close()
    args = "--matrix {} -c {} --outFileName {} -t {} --maximalDistance {} --distanceShortRange {} \
            --distanceLongRange {} --orderBy {} --dpi {} --legend"\
                .format(ROOT + 'test_matrix.scool',
                        ROOT + 'scHicPlotClusterProfiles/cluster_kmeans.txt', outfile.name,
                        1, 9000000, 2000000, 12000000,
                        'orderByFile', 300
                        ).split()
    scHicPlotClusterProfiles.main(args)
    test_image_path = ROOT + 'scHicPlotClusterProfiles/plot_maxDistance.png'
    res = compare_images(test_image_path, outfile.name, tolerance)
    assert res is None, res

    # os.unlink(outfile.name)


def test_version():
    args = "--version".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicPlotClusterProfiles.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0


def test_help():
    args = "--help".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicPlotClusterProfiles.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0
