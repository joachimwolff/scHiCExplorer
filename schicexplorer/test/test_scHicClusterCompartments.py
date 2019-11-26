import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp

from schicexplorer import scHicClusterCompartments

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
        --outFileName {} -t {}".format(ROOT + 'test_matrix.mcool',
                                       3, 'kmeans', outfile.name, 2).split()
    scHicClusterCompartments.main(args)
    assert are_files_equal(ROOT + "scHicClusterCompartments/cluster_kmeans.txt", outfile.name)


def test_spectral():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} ".format(ROOT + 'test_matrix.mcool',
                                        3, 'spectral', outfile.name, 2).split()
    scHicClusterCompartments.main(args)
    assert are_files_equal(ROOT + "scHicClusterCompartments/cluster_spectral.txt", outfile.name)


def test_kmeans_binarization():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} --binarization".format(ROOT + 'test_matrix.mcool',
                                                      3, 'kmeans', outfile.name, 2).split()
    scHicClusterCompartments.main(args)
    assert are_files_equal(ROOT + "scHicClusterCompartments/cluster_kmeans_binarization.txt", outfile.name)


def test_kmeans_histonmark():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} --histonMarkType {} --binarization --norm".format(ROOT + 'test_matrix.mcool',
                                                                                 3, 'kmeans', outfile.name, 2, ROOT + 'scHicClusterCompartments/mm9_H3K36me3.bed.gz').split()
    scHicClusterCompartments.main(args)
    assert are_files_equal(ROOT + "scHicClusterCompartments/cluster_kmeans_binarization_norm_histon_track.txt", outfile.name)


def test_spectral_histonmark():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} --histonMarkType {} --binarization --norm".format(ROOT + 'test_matrix.mcool',
                                                                                 3, 'spectral', outfile.name, 2, ROOT + 'scHicClusterCompartments/mm9_H3K36me3.bed.gz').split()
    scHicClusterCompartments.main(args)
    assert are_files_equal(ROOT + "scHicClusterCompartments/cluster_spectral_binarization_norm_histon_track.txt", outfile.name)


def test_spectral_extraTrack():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} --binarization --extraTrack {} --norm".format(ROOT + 'test_matrix.mcool',
                                                                             3, 'spectral', outfile.name, 2, ROOT + 'scHicClusterCompartments/mm9_gene.bed.gz').split()
    scHicClusterCompartments.main(args)
    assert are_files_equal(ROOT + "scHicClusterCompartments/cluster_spectral_binarization_norm_gene_track.txt", outfile.name)


def test_kmeans_extraTrack():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} --binarization --extraTrack {} --norm".format(ROOT + 'test_matrix.mcool',
                                                                             3, 'kmeans', outfile.name, 2, ROOT + 'scHicClusterCompartments/mm9_gene.bed.gz').split()
    scHicClusterCompartments.main(args)
    assert are_files_equal(ROOT + "scHicClusterCompartments/cluster_kmeans_binarization_norm_gene_track.txt", outfile.name)


def test_kmeans_norm():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} --binarization --norm".format(ROOT + 'test_matrix.mcool',
                                                             3, 'kmeans', outfile.name, 2).split()
    scHicClusterCompartments.main(args)
    assert are_files_equal(ROOT + "scHicClusterCompartments/cluster_kmeans_binarization_norm.txt", outfile.name)


def test_version():
    args = "--version".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicClusterCompartments.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0


def test_help():
    args = "--help".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicClusterCompartments.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0
