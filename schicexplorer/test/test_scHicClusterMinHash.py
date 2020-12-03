import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp
from matplotlib.testing.compare import compare_images

from schicexplorer import scHicClusterMinHash

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test-data/")
tolerance = 60


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
        --outFileName {} -t {} -nh {} -dim_pca 100".format(ROOT + 'test_matrix.scool',
                                                           3, 'kmeans', outfile.name, 3, 800).split()
    scHicClusterMinHash.main(args)
    assert are_files_equal_clustering(ROOT + "scHicClusterMinHash/cluster_kmeans.txt", outfile.name)


def test_agglomerative_ward():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -nh {} --noPCA".format(ROOT + 'test_matrix.scool',
                                                      3, 'agglomerative_ward', outfile.name, 3, 800).split()
    scHicClusterMinHash.main(args)
    assert are_files_equal_clustering(ROOT + "scHicClusterMinHash/cluster_kmeans.txt", outfile.name)


def test_agglomerative_complete():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -nh {} --noPCA ".format(ROOT + 'test_matrix.scool',
                                                       3, 'agglomerative_complete', outfile.name, 3, 800).split()
    scHicClusterMinHash.main(args)
    assert are_files_equal_clustering(ROOT + "scHicClusterMinHash/cluster_kmeans.txt", outfile.name)


def test_agglomerative_average():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -nh {} --noPCA ".format(ROOT + 'test_matrix.scool',
                                                       3, 'agglomerative_average', outfile.name, 3, 800).split()
    scHicClusterMinHash.main(args)
    assert are_files_equal_clustering(ROOT + "scHicClusterMinHash/cluster_kmeans.txt", outfile.name)


def test_agglomerative_single():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -nh {} --noPCA ".format(ROOT + 'test_matrix.scool',
                                                       3, 'agglomerative_single', outfile.name, 3, 800).split()
    scHicClusterMinHash.main(args)
    assert are_files_equal_clustering(ROOT + "scHicClusterMinHash/cluster_kmeans.txt", outfile.name)


def test_birch():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -nh {} -dim_pca 50".format(ROOT + 'test_matrix.scool',
                                                          3, 'birch', outfile.name, 3, 800).split()
    scHicClusterMinHash.main(args)
    assert are_files_equal_clustering(ROOT + "scHicClusterMinHash/cluster_kmeans.txt", outfile.name)


def test_spectral():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -nh {}".format(ROOT + 'test_matrix.scool',
                                              3, 'spectral', outfile.name, 2, 800).split()
    scHicClusterMinHash.main(args)
    assert are_files_equal_clustering(ROOT + "scHicClusterMinHash/cluster_spectral.txt", outfile.name)


def test_spectral_chromosomes():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -nh {} --chromosomes {} ".format(ROOT + 'test_matrix.scool',
                                                                3, 'spectral', outfile.name, 2, 800, "chr1 chr2").split()
    scHicClusterMinHash.main(args)
    assert are_files_equal_clustering(ROOT + "scHicClusterMinHash/cluster_spectral_chromosomes.txt", outfile.name)


# some issue with the test data, real world data works fine
@pytest.mark.xfail
def test_kmeans_euclidean():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -nh {} --euclideanModeMinHash ".format(ROOT + 'test_matrix.scool',
                                                                      3, 'kmeans', outfile.name, 2, 800).split()
    scHicClusterMinHash.main(args)
    assert are_files_equal_clustering(ROOT + "scHicClusterMinHash/cluster_kmeans_exact.txt", outfile.name)


# some issue with the test data, real world data works fine
# @pytest.mark.xfail
def test_spectral_euclidean():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile_plot = NamedTemporaryFile(prefix='pca_plot_', delete=False)
    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -nh {} --euclideanModeMinHash -csp {} --colorMap {} --dpi {} --fontsize {} --figuresize {} {}".format(ROOT + 'test_matrix.scool',
                                                                                                                                     3, 'spectral', outfile.name, 2, 800, outfile_plot.name,
                                                                                                                                     'tab10', 100, 5, 10, 5).split()
    scHicClusterMinHash.main(args)
    assert are_files_equal_clustering(ROOT + "scHicClusterMinHash/cluster_spectral_euclidean.txt", outfile.name)

    res = compare_images(ROOT + "scHicClusterMinHash/plot_pc1_pc2.png", outfile_plot.name + '.png', tolerance)
    assert res is None, res
    res = compare_images(ROOT + "scHicClusterMinHash/plot_pc2_pc3.png", outfile_plot.name + '.png', tolerance)
    assert res is None, res


def test_kmeans_saveMemory():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -nh {} -dim_pca 100 --saveMemory".format(ROOT + 'test_matrix.scool',
                                                                        3, 'kmeans', outfile.name, 3, 800).split()
    scHicClusterMinHash.main(args)
    assert are_files_equal_clustering(ROOT + "scHicClusterMinHash/cluster_kmeans.txt", outfile.name)


def test_agglomerative_single_saveMemory():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -nh {} --noPCA --saveMemory".format(ROOT + 'test_matrix.scool',
                                                                   3, 'agglomerative_single', outfile.name, 3, 800).split()
    scHicClusterMinHash.main(args)
    assert are_files_equal_clustering(ROOT + "scHicClusterMinHash/cluster_kmeans.txt", outfile.name)


def test_kmeans_euclidean_saveMemory():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -nh {} --euclideanModeMinHash --saveMemory ".format(ROOT + 'test_matrix.scool',
                                                                                   3, 'kmeans', outfile.name, 2, 800).split()
    scHicClusterMinHash.main(args)
    assert are_files_equal_clustering(ROOT + "scHicClusterMinHash/cluster_kmeans_exact.txt", outfile.name)


def test_version():
    args = "--version".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicClusterMinHash.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0


def test_help():
    args = "--help".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicClusterMinHash.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0
