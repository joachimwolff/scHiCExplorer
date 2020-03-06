scHiCExplorer
=============

Set of programs to process, normalize, analyse and visualize single-cell Hi-C data
----------------------------------------------------------------------------------


Availability
------------

scHiCExplorer is available as a **command line suite of tools** on this `GitHub repository <https://github.com/joachimwolff/scHiCExplorer>`_. scHiCExplorer is a general use single-cell Hi-C analysis software, to process raw single-cell Hi-C data we provide a demultiplexing tool for data provided by Nagano 2017. For all other protocols the demultiplexing must be computed by a third party tool. However, as long as per cell one forward and reverse FASTQ respectivly after mapping a BAM/SAM file is provided, scHiCExplorer is able to process it.


The following is the list of tools available in scHiCExplorer
-------------------------------------------------------------


=================================== ===============================================================================
tool                                description
=================================== ===============================================================================
:ref:`scHicDemultiplex`             Demultiplexes the samples by their barcodes to one FASTQ file per samples
:ref:`scHicMergeToMCool`            Merges all single-cell Hi-C matrices to one
:ref:`scHicMergeMatrixBins`         Changes the resolution of the matrices 
:ref:`scHicQualityControl`          Estimates the quality of scHi-C datasets
:ref:`scHicAdjustMatrix`            Keeps / removes chromosomes / contigs / scaffolds of all samples 
:ref:`scHicNormalize`               Normalizes the read coverage of all samples to the lowest read coverage
:ref:`scHicCorrectMatrices`         Corrects all samples with Knight-Ruiz correction
:ref:`scHicInfo`                    Retrieve information about the mcool matrix
:ref:`scHicCluster`                 Cluster all samples on raw data or uses dimension reduction knn or pca 
:ref:`scHicClusterMinHash`          Cluster all samples on knn computed by approximate knn search via LSH
:ref:`scHicClusterSVL`              Cluster all samples based on short vs long range contact ratio  
:ref:`scHicClusterCompartments`     Cluster all samples based on A / B scHicClusterCompartments  
:ref:`scHicConsensusMatrices`       Computes the consensus matrices based on clustering
:ref:`scHicPlotClusterProfiles`     Plots the cluster profiles with all samples
:ref:`scHicPlotConsensusMatrices`   Estimates the quality of Hi-C dataset
=================================== ===============================================================================

Getting Help
------------

* For all kind of questions, suggesting changes/enhancements and to report bugs, please create an issue on `our GitHub repository <https://github.com/joachimwolff/scHiCExplorer>`_

Contents:
---------

.. toctree::
   :maxdepth: 2

   content/installation
   content/list-of-tools
   content/News
   content/Example_analysis



This tool suite is developed by Joachim Wolff from the `Bioinformatics Lab <http://bioinf.uni-freiburg.de/>`_ of the `Albert-Ludwigs-University Freiburg <http://www.uni-freiburg.de>`_, Germany.
