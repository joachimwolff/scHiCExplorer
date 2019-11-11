.. _Example_analysis:

Analysis of single-cell Hi-C data
=================================

The analysis of single-cell Hi-C data deals is partially similar to regular Hi-C data analysis, the pre-processing of data i.e. mapping and the creation
of a Hi-C interaction matrix and the correction of the data works can be adapted from a Hi-C data anlysis. However, single-cell Hi-C deals
with different issues as the be mentioned the pre-processing of the fastq files (demultiplexing) to assoziate the reads to one sample (to one cell). 
Everything that applies to Hi-C also applies to single-cell Hi-C exepct in single-cell the work is done with a few thousand cells and not one. Additonal, the read coverage
in single-cell Hi-C is not in the hundereds of millions but on the used data form Nagano 2017 on average 1.5 million. This leads to other necessary treatments in the quality 
control of the samples and a need for a normalization to a equal read coverage for all cells.


In this tutorial we work with the 'diploid' data from Nagano 2017 (GSE94489). 

**Disclaimer**

The raw fastq data is around 1,04 TB of size and the download speed is limited to a few Mb/s by NCBI. To decrease the download time it is recommended to download the files in parallel if enough disk space is available.
Furthermore, please consider the data needs to be demultiplexed and mapped which needs additional disk space.

If you do not want to download, demultiplex, map and build the matrices on your own, a precomputed mcool matrix is provided #TODO: here.

Download of the fastq files
^^^^^^^^^^^^^^^^^^^^^^^^^^^

As the very first step the raw, non-demultiplexed fastq files need to be downloaded. Please downloade the files directly from NCBI GEO and not e.g. from EMBL ENA, we 
have seen that these files miss the barcode information.

To download the fastq files the SRR sample number must be known, for not all samples only one SRR number was given, these samples were therefore not included in this tutorial.

.. code-block:: bash

    SRR5229019	GSM2476401	Diploid_11
    SRR5229021	GSM2476402	Diploid_12
    SRR5229023	GSM2476403	Diploid_13
    SRR5229025	GSM2476404	Diploid_15_16_17_18
    SRR5229027	GSM2476405	Diploid_1_6
    SRR5229029	GSM2476406	Diploid_19
    SRR5229031	GSM2476407	Diploid_20
    SRR5229033	GSM2476408	Diploid_2_14
    SRR5229035	GSM2476409	Diploid_21
    SRR5229037	GSM2476410	Diploid_22
    SRR5229039	GSM2476411	Diploid_23
    SRR5229041	GSM2476412	Diploid_24
    SRR5229043	GSM2476413	Diploid_25
    SRR5229045	GSM2476414	Diploid_26
    SRR5229047	GSM2476415	Diploid_3
    SRR5229049	GSM2476416	Diploid_4
    SRR5229051	GSM2476417	Diploid_5_10
    SRR5229053	GSM2476418	Diploid_7
    SRR5229055	GSM2476419	Diploid_8
    SRR5229057	GSM2476420	Diploid_9
    SRR5507552	GSM2598387	Diploid_28_29
    SRR5507553	GSM2598388	Diploid_30_31
    SRR5507554	GSM2598389	Diploid_32_33
    SRR5507555	GSM2598390	Diploid_34_35

Excluded: GSM2598386 / Diploid_27


Download the each file via:

.. code-block:: bash

    $ fastq-dump SRR5229019

Alternatively, download all with one command:

.. code-block:: bash

    $ echo SRR5229019,SRR5229021,SRR5229023,SRR5229025,SRR5229027,SRR5229031,SRR5229033,SRR5229035,SRR5229037,SRR5229039,SRR5229041,SRR5229043,SRR5229045,SRR5229047,SRR5229049,SRR5229051,SRR5229053,SRR5229055,SRR5229057,SRR5507553,SRR5507554,SRR5507555 |  sed "s/,/\n/g" | xargs -n1 -P 22 -I {} sh -c "fastq-dump {}" 



Demultiplexing
^^^^^^^^^^^^^^

Each downloaded file needs to be demultiplexed. To do so the barcodes per sample (https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE94489&format=file&file=GSE94489%5FREADME%2Etxt) and the SRR to sample mapping needs to be provided:


.. code-block:: bash

    $ scHicDemultiplex -f "FASTQ_FILE" --srrToSampleFile samples.txt --barcodeFile GSE94489_README.txt --threads 20


scHicDemultiplex creates a folder 'demultiplexed' containing the demultiplexed fastq files splitted as forward and reverse reads and follows the scheme:

.. code-block::

    sample_id_barcode_RX.fastq.gz

For example:

.. code-block::

    Diploid_15_AGGCAGAA_CTCTCTAT_R1.fastq.gz


Please consider that the time to demultiplex the file SRR5229025, which itself is 4.1 GB takes already ~35 mins, to demultiplex the full 1 TB dataset will take around 6 days to compute.


Mapping
^^^^^^^

After demultiplexing, each forward and reverse strand file needs to be mapped as usual in Hi-C as single-paired files. Foe this tutorial we use bowtie2 and the mm10 index:

.. code-block:: bash

    $ wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm10.zip
    $ mkdir mm10 && unzip mm10.zip -d mm10


.. code-block:: bash

    $ bowtie2 -x mm10/mm10 --threads 8 -U ../original_data/SRR1956527_1.fastq.gz --reorder | samtools view -Shb - > SRR1956527_1.bam
    $ ls demultiplexed |  xargs -n1 -P 5 -I {} sh -c "bowtie2 -x mm10/mm10 --threads 5 -U demultiplexed/{} --reorder | samtools view -Shb - > {}.bam"



Creation of Hi-C interaction matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As a last step, the matrices for each cell need to be created, we use the tool 'hicBuildMatrix' from HiCExplorer:

.. code-block:: bash

    $ ls *.bam |  tr '\n' ' ' | xargs -n 2 -P 1 -d ' ' | xargs -n1 -P5 -I {} bash -c 'multinames=$1;outname=$(echo $multinames | cut -d" " -f 1 | sed -r "s?(^.*)_R[12]\..*?\\1?"); mkdir ${outname}_QC && hicBuildMatrix -s $multinames --binSize 1000000 --QCfolder  ${outname}_QC -o ${outname}.cool --threads 4' -- {}


To make this step more automated, it is recommend to use either a platform like hicexplorer.usegalaxy.eu or to use a batch script:

.. code-block:: bash
    
    $ ls -1 *.bam |  xargs -n2 -P 1 -I {} sh -c "hicBuildMatrix -s {} {} --binSize 1000000 name.cool"


After the Hi-C interaction matrices for each cell is created, the matrices are pooled together to one mcool matrix:

.. code-block:: bash

    $ scHicMergeMatrixBins --matrices matrices/* --outFileName nagano2017_raw.mcool


Quality control
^^^^^^^^^^^^^^^

Quality control is the crucuial step in preprocessing of all HTS related data. For single-cell experiements the read covarage 
per sample needs to be on a minmal level, and all matrices needs to be not broken and contain all the same chromosomes. Especially the last two issues are 
likely to rise in single-cell Hi-C data because the read coverage is with around 1 million reads, in contrast to regular Hi-C with a few 
hundert million, quite low and therefore it is more likley that simply no data for small chromosomes is present. 
To guarantee these requirements the quality control works in three steps: 

1. Only matrices which contain all listed chromosomes are accepted
2. Only matrices which have a minimum read coverage are accepted
3. The matrix must have a minum denisity of recorded data points close to the main diagonal.

.. code-block:: bash

    $ scHicQualityControl --matrix nagano2017_raw.mcool --outputMcool nagano2017_qc.mcool --minimumReadCoverage 1000000 --minimumDensity 0.001 --maximumRegionToConsider 30000000 --outFileNameReadCoverage read_coverage.png --outFileNameSparsity sparsity.png --chromosomes chr1 chr2 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19

For this tutorial a minimum read coverage of 1 million and a density of 0.1% is used in range of 30MB around the main diagonal. The above command creates certain files:

1. A mcool matrix containing only samples with matrices that passed the quality settings.
2. A plot showing the sparsity of all samples. Use this plot to adjust the minimumDensity parameter.
3. A plot showing the read coverage of all samples, use this plot to adjust the minimum read coverage parameter.
4. Three files containing information about which samples were removed / kept.
5. An HTML report given additional quality control information.


Analysis
^^^^^^^^

The analysis of single-cell Hi-C data investigates the chromatin folding changes during the cell cycle. 
To compute this, the clustering of the cells and a correct ordering within a cluster is the key step for this analysis.

scHiCExplorer uses a flatting approach to create out of the two dimensional 2D interaction matrices a one dimensional vector to have in the end 
a numer of samples times number of bins^2 matrix. For example: Nagano 2017 has around 3000 cells and using a 1MB binning approach results for the mouse genome in
2600 times 2600 matrix. After flattening, the matrix which is used to operate on is 3000 * (2600 * 2600) = 3000 * 6760000. 

Two aproaches to apply clustering are now possible: 

1. Compute the clustering directly on the matrix.
2. Reduce the dimensions first and apply clustering.

Option one works if the resolution of the interaction matrices are not too high, i.e. 1MB leads to 6.7 million features which is already a lot, but todays computers can handle this.
However, it looks different if the resolution is increased to e.g. regular Hi-C matrix resolution of 10kb. In this case the binned matrix is not 2600 * 2600, but 260000 * 260000 which is 67.6 billion.
To work on such many features would be problematic in terms of computational time and, it is questionable if a computer with enough main memory is available.
To overcome this, a dimension reduction is necessary. To reduce the number of dimensions scHiCExplorer provides three approaches: MinHash, SVL and Compartments.

The first approach uses a local sensitive hashing approach to compute the nearest neighbors, with it, it reduces the number of dimensions to the number of samples where each entry represents how close the samples are. 
Approach two, SVL for short vs long distances, computes per chromosome the ratio of the sum of short range contancts vs. the sum of long range contacts, the number of dimensions is therefore reduced to the number of to be considered chromosomes. 
Approach number three, compartments, computes the A/B compartments per chromosome and reduces the number of dimensions to the square root.

In the following, all four approaches are shown.

