import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp

from schicexplorer import scHicDemultiplex
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test-data/")


def test_demultiplexing():
    output_folder = mkdtemp(prefix="output_demultiplex")

    args = "--fastq {} --barcodeFile {} --srrToSampleFile {} \
        --outputFolder {} -t {} --bufferSize {}".format(ROOT + 'scHicDemultiplex/SRR5229025.fastq.gz',
                                                        ROOT + 'scHicDemultiplex/GSE94489_README.txt',
                                                        ROOT + 'scHicDemultiplex/samples.txt',
                                                        output_folder,
                                                        1,
                                                        1000).split()
    scHicDemultiplex.main(args)

    assert set(os.listdir(ROOT + '/scHicDemultiplex/demultiplexed/')) == set(os.listdir(output_folder))


def test_version():
    args = "--version".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicDemultiplex.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0


def test_help():
    args = "--help".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicDemultiplex.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0
