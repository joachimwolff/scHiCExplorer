import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp

import cooler
from schicexplorer import scHicCorrectMatrices
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test-data/")


# def are_files_equal(file1, file2, delta=2, skip=0):
#     equal = True
#     if delta:
#         mismatches = 0
#     with open(file1) as textfile1, open(file2) as textfile2:
#         for i, (x, y) in enumerate(zip(textfile1, textfile2)):
#             if i < skip:
#                 continue
#             if x != y:
#                 if delta:
#                     mismatches += 1
#                     if mismatches > delta:
#                         equal = False
#                         break
#                 else:
#                     equal = False
#                     break
#     return equal

def test_kmeans():
    outfile = NamedTemporaryFile(suffix='.mcool', delete=False)

    outfile.close()
    args = "--matrix {} --outFileName {} -t {} ".format(ROOT + 'test_matrix.mcool',
                                outfile.name).split()
    scHicCorrectMatrices.main(args)


    test = hm.hiCMatrix(
        ROOT + "hicCorrectMatrix/small_test_matrix_KRcorrected_chrUextra_chr3LHet.h5")
    new = hm.hiCMatrix(outfile.name)
    nt.assert_almost_equal(test.matrix.data, new.matrix.data, decimal=5)
    nt.assert_equal(test.cut_intervals, new.cut_intervals)

    os.unlink(outfile.name)


    assert 

