import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp
import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from schicexplorer import scHicPlotConsensusMatrices
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test-data/")

tolerance = 30


def test_plot():
    outfile = NamedTemporaryFile(suffix='.png', delete=False)

    outfile.close()
    args = "--matrix {} --outFileName {} -t {} --dpi {}".format(ROOT + 'scHicConsensusMatrices/consensus_matrix.mcool',
                                                                outfile.name, 1, 300
                                                                ).split()
    scHicPlotConsensusMatrices.main(args)
    test_image_path = ROOT + 'scHicPlotConsensusMatrices/plot.png'
    res = compare_images(test_image_path, outfile.name, tolerance)
    assert res is None, res

    os.unlink(outfile.name)


def test_plot_chr():
    outfile = NamedTemporaryFile(suffix='.png', delete=False)

    outfile.close()
    args = "--matrix {} --outFileName {} -t {} --dpi {} -c {}".format(ROOT + 'scHicConsensusMatrices/consensus_matrix.mcool',
                                                                      outfile.name, 1, 300, "chr1"
                                                                      ).split()
    scHicPlotConsensusMatrices.main(args)
    test_image_path = ROOT + 'scHicPlotConsensusMatrices/plot_chr1.png'
    res = compare_images(test_image_path, outfile.name, tolerance)
    assert res is None, res

    os.unlink(outfile.name)


def test_plot_multi_chr():
    outfile = NamedTemporaryFile(suffix='.png', delete=False)

    outfile.close()
    args = "--matrix {} --outFileName {} -t {} --dpi {} -c {}".format(ROOT + 'scHicConsensusMatrices/consensus_matrix.mcool',
                                                                      outfile.name, 1, 300, 'chr1 chr2'
                                                                      ).split()
    scHicPlotConsensusMatrices.main(args)
    test_image_path = ROOT + 'scHicPlotConsensusMatrices/plot_chr1_chr2.png'
    res = compare_images(test_image_path, outfile.name, tolerance)
    assert res is None, res

    os.unlink(outfile.name)


def test_version():
    args = "--version".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicPlotConsensusMatrices.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0


def test_help():
    args = "--help".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicPlotConsensusMatrices.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0
