import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp

from schicexplorer import scHicClusterMinHash

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test-data/")


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

def test_kmeans():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -h {}".format(ROOT + 'test_matrix.mcool',
                                3, 'kmeans', outfile.name, 1, 800).split()
    scHicClusterMinHash.main(args)
    assert are_files_equal(ROOT + "scHicClusterMinHash/cluster_kmeans.txt", outfile.name)

def test_spectral():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -h {}".format(ROOT + 'test_matrix.mcool',
                                3, 'spectral', outfile.name, 1, 800).split()
    scHicClusterMinHash.main(args)
    assert are_files_equal(ROOT + "scHicClusterMinHash/cluster_spectral.txt", outfile.name)

def test_spectral_chromosomes():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -h {} --chromosomes {} ".format(ROOT + 'test_matrix.mcool',
                                3, 'spectral', outfile.name, 1, 800, "chr1 chr2").split()
    scHicClusterMinHash.main(args)
    assert are_files_equal(ROOT + "scHicClusterMinHash/cluster_spectral_chromosomes.txt", outfile.name)

def test_kmeans_exact():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -h {}".format(ROOT + 'test_matrix.mcool',
                                3, 'kmeans', outfile.name, 1, 800).split()
    scHicClusterMinHash.main(args)
    assert are_files_equal(ROOT + "scHicClusterMinHash/cluster_kmeans_exact.txt", outfile.name)


def test_spectral_exact():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -h {}".format(ROOT + 'test_matrix.mcool',
                                3, 'spectral', outfile.name, 1, 800).split()
    scHicClusterMinHash.main(args)
    assert are_files_equal(ROOT + "scHicClusterMinHash/cluster_spectral_exact.txt", outfile.name)