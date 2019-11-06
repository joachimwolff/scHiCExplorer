scHiCExplorer tools
===================

.. contents::
    :local:


+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
| tool                                 | type             | input files                       | main output file(s)                         | application                                                                       |
+======================================+==================+===================================+=============================================+===================================================================================+
|:ref:`scHicAdjustMatrix`              | preprocessing    | 1 genome FASTA file               | bed file with restriction site coordinates  | Identifies the genomic locations of restriction sites                             |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicAdjustMatrix`              | preprocessing    | 1 genome FASTA file               | bed file with restriction site coordinates  | Identifies the genomic locations of restriction sites                             |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicAdjustMatrix`              | preprocessing    | 1 genome FASTA file               | bed file with restriction site coordinates  | Identifies the genomic locations of restriction sites                             |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicAdjustMatrix`              | preprocessing    | 1 genome FASTA file               | bed file with restriction site coordinates  | Identifies the genomic locations of restriction sites                             |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicAdjustMatrix`              | preprocessing    | 1 genome FASTA file               | bed file with restriction site coordinates  | Identifies the genomic locations of restriction sites                             |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicAdjustMatrix`              | preprocessing    | 1 genome FASTA file               | bed file with restriction site coordinates  | Identifies the genomic locations of restriction sites                             |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicAdjustMatrix`              | preprocessing    | 1 genome FASTA file               | bed file with restriction site coordinates  | Identifies the genomic locations of restriction sites                             |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicAdjustMatrix`              | preprocessing    | 1 genome FASTA file               | bed file with restriction site coordinates  | Identifies the genomic locations of restriction sites                             |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicAdjustMatrix`              | preprocessing    | 1 genome FASTA file               | bed file with restriction site coordinates  | Identifies the genomic locations of restriction sites                             |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicAdjustMatrix`              | preprocessing    | 1 genome FASTA file               | bed file with restriction site coordinates  | Identifies the genomic locations of restriction sites                             |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicAdjustMatrix`              | preprocessing    | 1 genome FASTA file               | bed file with restriction site coordinates  | Identifies the genomic locations of restriction sites                             |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicAdjustMatrix`              | preprocessing    | 1 genome FASTA file               | bed file with restriction site coordinates  | Identifies the genomic locations of restriction sites                             |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicAdjustMatrix`              | preprocessing    | 1 genome FASTA file               | bed file with restriction site coordinates  | Identifies the genomic locations of restriction sites                             |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`scHicAdjustMatrix`              | preprocessing    | 1 genome FASTA file               | bed file with restriction site coordinates  | Identifies the genomic locations of restriction sites                             |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+


General principles
^^^^^^^^^^^^^^^^^^

A typical HiCExplorer command could look like this:

.. code:: bash

 $ hicPlotMatrix -m myHiCmatrix.h5 \
 -o myHiCmatrix.pdf \
 --clearMaskedBins \
 --region chrX:10,000,000-15,000,000 \
 --vMin -4 --vMax 4 \


You can always see all available command-line options via --help:

.. code:: bash

 $ hicPlotMatrix --help


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