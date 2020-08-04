import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp

from schicexplorer import scHicClusterSVL
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test-data/")


def are_files_equal_clustering(file1, file2, number_of_clusters=3, delta=2, skip=0):
    equal = True
    if delta:
        mismatches = 0
    numberOfClusters = set()
    with open(file1) as textfile1, open(file2) as textfile2:
        for i, (x, y) in enumerate(zip(textfile1, textfile2)):
            if i < skip:
                continue
            x = x.split(' ')
            y = y.split(' ')
            numberOfClusters.add(y[1])
            x[0] = x[0].lstrip('/cells/')
            y[0] = y[0].lstrip('/cells/')

            if x[0] != y[0]:
                if delta:
                    mismatches += 1
                    if mismatches > delta:
                        equal = False
                        break
                else:
                    equal = False
                    break
    if len(numberOfClusters) == number_of_clusters:
        return equal
    else:
        return False
    return equal


def test_kmeans():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -ds  {} -dl {} ".format(ROOT + 'test_matrix.scool',
                                                       3, 'kmeans', outfile.name, 2, 2000000, 12000000).split()
    scHicClusterSVL.main(args)
    assert are_files_equal_clustering(ROOT + "scHicClusterSVL/cluster_kmeans.txt", outfile.name)


def test_spectral():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -ds  {} -dl {} ".format(ROOT + 'test_matrix.scool',
                                                       3, 'spectral', outfile.name, 2, 2000000, 12000000).split()
    scHicClusterSVL.main(args)
    assert are_files_equal_clustering(ROOT + "scHicClusterSVL/cluster_spectral.txt", outfile.name)


def test_version():
    args = "--version".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicClusterSVL.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0


def test_help():
    args = "--help".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicClusterSVL.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0
