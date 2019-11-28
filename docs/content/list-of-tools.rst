scHiCExplorer tools
===================

.. contents::
    :local:


+--------------------------------------+------------------+----------------------------------------+----------------------------------------------+-----------------------------------------------------------------------------------+
| tool                                 | type             | input files                            | main output file(s)                          | application                                                                       |
+======================================+==================+========================================+==============================================+===================================================================================+
|:ref:`scHicDemultiplex`               | preprocessing    | 1 FASTQ file                           | 2*n demultiplexed FASTQ files                | Demultiplexes the samples by their barcodes to one FASTQ file per samples         |
|                                      |                  | SRR to sample mapping file             |                                              |                                                                                   |
|                                      |                  | Barcode file                           |                                              |                                                                                   |
+--------------------------------------+------------------+----------------------------------------+----------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicMergeToMCool`              | preprocessing    | n Hi-C matrices in cool format         | One mcool file containg all Hi-C matrices    | Merges all single-cell Hi-C matrices to one                                       |
+--------------------------------------+------------------+----------------------------------------+----------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicMergeMatrixBins`           | preprocessing    | mcool Hi-C matrix                      | mcool Hi-C matrix                            | Changes the resolution of the matrices                                            |
+--------------------------------------+------------------+----------------------------------------+----------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicQualityControl`            | preprocessing    | mcool Hi-C matrix                      | One mcool file, two qc images, qc report     | Checks the quality of all samples and removes bad ones                            |
+--------------------------------------+------------------+----------------------------------------+----------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicAdjustMatrix`              | preprocessing    | mcool Hi-C matrix                      | mcool Hi-C matrix                            | Keeps / removes chromosomes / contigs / scaffolds of all samples                  |
+--------------------------------------+------------------+----------------------------------------+----------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicNormalize`                 | preprocessing    | mcool Hi-C matrix                      | mcool Hi-C matrix                            | Normalizes the read coverage of all samples to the lowest read coverage           |
+--------------------------------------+------------------+----------------------------------------+----------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicCorrectMatrices`           | preprocessing    | mcool Hi-C matrix                      | mcool Hi-C matrix                            | Corrects all samples with Knight-Ruiz correction                                  |
+--------------------------------------+------------------+----------------------------------------+----------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicInfo   `                   | information      | mcool Hi-C matrix                      | information about the mcool matrix           | Retrieve information about the mcool matrix: resolution, number of samples, etc   |
+--------------------------------------+------------------+----------------------------------------+----------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicCluster`                   | analysis         | mcool Hi-C matrix                      | text file with sample to cluster association | Cluster all samples on raw data or uses dimension reduction knn or pca            |
+--------------------------------------+------------------+----------------------------------------+----------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicClusterMinHash`            | analysis         | mcool Hi-C matrix                      | text file with sample to cluster association | Cluster all samples on knn computed by approximate knn search via LSH             |
+--------------------------------------+------------------+----------------------------------------+----------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicClusterSVL`                | analysis         | mcool Hi-C matrix                      | text file with sample to cluster association | Cluster all samples based on short vs long range contact ratio                    |
+--------------------------------------+------------------+----------------------------------------+----------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicClusterCompartments`       | analysis         | mcool Hi-C matrix                      | text file with sample to cluster association | Cluster all samples based on A / B scHicClusterCompartments                       |
|                                      |                  | (gene or histon track)                 |                                              |                                                                                   | 
+--------------------------------------+------------------+----------------------------------------+----------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicConsensusMatrix`           | analysis         | mcool Hi-C matrix,                     | mcool Hi-C matrix with consensus matrices    | Computes the consensus matrices based on clustering                               |
|                                      |                  | txt file sample to cluster association |                                              |                                                                                   |
+--------------------------------------+------------------+----------------------------------------+----------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicPlotClusterProfiles`       | visualization    | mcool Hi-C matrix                      | one image with cluster profiles              | Plots the cluster profiles with all samples                                       |
|                                      |                  | txt file sample to cluster association |                                              |                                                                                   |
+--------------------------------------+------------------+----------------------------------------+----------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicPlotConsensusMatrices`     | visualization    | mcool Hi-C matrix                      | one image with consensus matrices            | Plots the cluster consensus matrices                                              |
|                                      |                  | txt file sample to cluster association |                                              |                                                                                   |
+--------------------------------------+------------------+----------------------------------------+----------------------------------------------+-----------------------------------------------------------------------------------+


General principles
^^^^^^^^^^^^^^^^^^

A typical scHiCExplorer command could look like this:

.. code:: bash

 $ scHicPlotClusterProfiles -m matrix.mcool \
 -o cluster_profiles.png \
 -c computed_clusters.txt \ 
 --dpi 300


You can always see all available command-line options via --help:

.. code:: bash

 $ scHicInfo -m matrix.mcool


Tools for demultiplexing of raw fastq files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:ref:`scHicDemultiplex`
"""""""""""""""""""""""


Tools for single-cell Hi-C data pre-preprocessing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:ref:`scHicAdjustMatrix`
""""""""""""""""""""""""
:ref:`scHicCorrectMatrices`
"""""""""""""""""""""""""""
:ref:`scHicMergeMatrixBins`
"""""""""""""""""""""""""""
:ref:`scHicMergeToMCool`
""""""""""""""""""""""""
:ref:`scHicNormalize`
"""""""""""""""""""""

Tools for information about the single-cell Hi-C matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:ref:`scHicInfo`
""""""""""""""""

Tools for single-cell Hi-C QC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:ref:`scHicQualityControl`
""""""""""""""""""""""""""

Tools for Hi-C data analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:ref:`scHicCluster`
"""""""""""""""""""
:ref:`scHicClusterCompartments`
"""""""""""""""""""""""""""""""
:ref:`scHicClusterMinHash`
""""""""""""""""""""""""""
:ref:`scHicClusterSVL`
""""""""""""""""""""""
:ref:`scHicConsensusMatrix`
"""""""""""""""""""""""""""

Tools for single-cell Hi-C visualization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:ref:`scHicPlotClusterProfiles`
"""""""""""""""""""""""""""""""
:ref:`scHicPlotConsensusMatrices`
"""""""""""""""""""""""""""""""""