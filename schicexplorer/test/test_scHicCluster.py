import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp
from matplotlib.testing.compare import compare_images

from schicexplorer import scHicCluster
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


# @pytest.mark.xfail
def test_kmeans():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {}".format(ROOT + 'test_matrix.scool',
                                       3, 'kmeans', outfile.name, 1).split()
    scHicCluster.main(args)
    assert are_files_equal_clustering(ROOT + "scHicCluster/cluster_kmeans.txt", outfile.name)


# @pytest.mark.xfail
def test_kmeans_chromosomes():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} --chromosomes {}".format(ROOT + 'test_matrix.scool',
                                                        3, 'kmeans', outfile.name, 1, "chr1 chr2").split()
    scHicCluster.main(args)
    assert are_files_equal_clustering(ROOT + "scHicCluster/cluster_kmeans_chromosomes.txt", outfile.name)


# @pytest.mark.xfail
def test_spectral():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {}".format(ROOT + 'test_matrix.scool',
                                       3, 'spectral', outfile.name, 1).split()
    scHicCluster.main(args)
    assert are_files_equal_clustering(ROOT + "scHicCluster/cluster_spectral.txt", outfile.name)


# @pytest.mark.xfail
def test_spectral_knn():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -drm {}".format(ROOT + 'test_matrix.scool',
                                               3, 'spectral', outfile.name, 1, "knn").split()
    scHicCluster.main(args)
    assert are_files_equal_clustering(ROOT + "scHicCluster/cluster_spectral_knn.txt", outfile.name)


# @pytest.mark.xfail
def test_spectral_pca():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -drm {}".format(ROOT + 'test_matrix.scool',
                                               3, 'spectral', outfile.name, 4, "pca").split()
    scHicCluster.main(args)
    assert are_files_equal_clustering(ROOT + "scHicCluster/cluster_spectral_pca.txt", outfile.name)


# @pytest.mark.xfail
def test_kmeans_knn():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -drm {}".format(ROOT + 'test_matrix.scool',
                                               3, 'kmeans', outfile.name, 2, "knn").split()
    scHicCluster.main(args)
    assert are_files_equal_clustering(ROOT + "scHicCluster/cluster_kmeans_knn.txt", outfile.name)


# @pytest.mark.xfail
def test_kmeans_pca():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -drm {}".format(ROOT + 'test_matrix.scool',
                                               3, 'kmeans', outfile.name, 3, "pca").split()
    scHicCluster.main(args)
    assert are_files_equal_clustering(ROOT + "scHicCluster/cluster_kmeans_pca.txt", outfile.name)


def test_kmeans_clustering():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {}".format(ROOT + 'test_matrix.scool',
                                       3, 'kmeans', outfile.name, 1).split()
    scHicCluster.main(args)
    assert are_files_equal_clustering(ROOT + "scHicCluster/cluster_kmeans.txt", outfile.name)


def test_kmeans_chromosomes_clustering():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} --chromosomes {}".format(ROOT + 'test_matrix.scool',
                                                        3, 'kmeans', outfile.name, 1, "chr1 chr2").split()
    scHicCluster.main(args)
    assert are_files_equal_clustering(ROOT + "scHicCluster/cluster_kmeans_chromosomes.txt", outfile.name)


def test_spectral_clustering_csp():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile_plot = NamedTemporaryFile(prefix='pca_plot_', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -csp {} --colorMap {} --dpi {} --fontsize {} --figuresize {} {}".format(ROOT + 'test_matrix.scool',
                                                                                                       3, 'spectral', outfile.name, 1, outfile_plot.name,
                                                                                                       'tab10', 100, 5, 15, 8).split()
    scHicCluster.main(args)
    assert are_files_equal_clustering(ROOT + "scHicCluster/cluster_spectral.txt", outfile.name)
    res = compare_images(ROOT + "scHicCluster/plot_pc1_pc2.png", outfile_plot.name + '.png', tolerance)
    assert res is None, res
    res = compare_images(ROOT + "scHicCluster/plot_pc2_pc3.png", outfile_plot.name + '.png', tolerance)
    assert res is None, res


# @pytest.mark.xfail
def test_spectral_knn_clustering():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -drm {}".format(ROOT + 'test_matrix.scool',
                                               3, 'spectral', outfile.name, 1, "knn").split()
    scHicCluster.main(args)
    assert are_files_equal_clustering(ROOT + "scHicCluster/cluster_spectral_knn.txt", outfile.name)


def test_spectral_pca_clustering_csp():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile_plot = NamedTemporaryFile(prefix='pca_plot_', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -drm {} -csp {} --colorMap {} --dpi {} --fontsize {} --figuresize {} {}".format(ROOT + 'test_matrix.scool',
                                                                                                               3, 'spectral', outfile.name, 4, "pca", outfile_plot.name,
                                                                                                               'tab10', 100, 5, 15, 8).split()
    scHicCluster.main(args)
    assert are_files_equal_clustering(ROOT + "scHicCluster/cluster_spectral_pca.txt", outfile.name)
    assert are_files_equal_clustering(ROOT + "scHicCluster/cluster_spectral.txt", outfile.name)
    res = compare_images(ROOT + "scHicCluster/plot_pc1_pc2.png", outfile_plot.name + '.png', tolerance)
    assert res is None, res
    res = compare_images(ROOT + "scHicCluster/plot_pc2_pc3.png", outfile_plot.name + '.png', tolerance)
    assert res is None, res


# @pytest.mark.xfail
def test_kmeans_knn_clustering_intra_chromosomal_plot():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile_plot = NamedTemporaryFile(prefix='pca_plot_', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -drm {} -ic -csp {} --colorMap {} --dpi {} --fontsize {} --figuresize {} {}".format(ROOT + 'test_matrix.scool',
                                                                                                                   3, 'kmeans', outfile.name, 2, "knn", outfile_plot.name,
                                                                                                                   'tab10', 100, 5, 15, 8).split()
    scHicCluster.main(args)
    assert are_files_equal_clustering(ROOT + "scHicCluster/cluster_kmeans_ic_knn.txt", outfile.name)

    res = compare_images(ROOT + "scHicCluster/plot_pc1_pc2.png", outfile_plot.name + '.png', tolerance)
    assert res is None, res
    res = compare_images(ROOT + "scHicCluster/plot_pc2_pc3.png", outfile_plot.name + '.png', tolerance)
    assert res is None, res


def test_kmeans_pca_clustering():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)

    outfile.close()
    args = "--matrix {} --numberOfClusters {} --clusterMethod {} \
        --outFileName {} -t {} -drm {} ".format(ROOT + 'test_matrix.scool',
                                                3, 'kmeans', outfile.name, 3, "pca").split()
    scHicCluster.main(args)
    assert are_files_equal_clustering(ROOT + "scHicCluster/cluster_kmeans_pca.txt", outfile.name)


def test_version():
    args = "--version".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicCluster.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0


def test_help():
    args = "--help".split()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        scHicCluster.main(args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0
