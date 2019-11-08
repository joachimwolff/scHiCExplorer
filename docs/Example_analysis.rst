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

To download the data can take a week or two, after demultiplexing and mapping, you need ~ 15 TB of free disk space.


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

.. code-block:: bash

    SRR5229059	GSM2476421	Haploid_10
    SRR5229061	GSM2476422	Haploid_11
    SRR5229063	GSM2476423	Haploid_12
    SRR5229065	GSM2476424	Haploid_13
    SRR5229067	GSM2476425	Haploid_14
    SRR5229084	GSM2476427	Haploid_16
    SRR5229086	GSM2476428	Haploid_17
    SRR5229088	GSM2476429	Haploid_18
    SRR5229090	GSM2476430	Haploid_19
    SRR5229103	GSM2476432	Haploid_20
    SRR5229105	GSM2476433	Haploid_2
    SRR5229107	GSM2476434	Haploid_3
    SRR5229111	GSM2476435	Haploid_4
    SRR5229113	GSM2476436	Haploid_5
    SRR5229115	GSM2476437	Haploid_6
    SRR5229117	GSM2476438	Haploid_7
    SRR5229119	GSM2476439	Haploid_8
    SRR5229121	GSM2476440	Haploid_9

Excluded: GSM2476426 / Haploid_15 and GSM2476431 / Haploid_1

Download the each file via:

.. code-block:: bash

    $ fastq-dump SRR5229019

Alternatively, download all with one command:

.. code-block:: bash

    $ echo SRR5229019,SRR5229021,SRR5229023,SRR5229025,SRR5229027,SRR5229031,SRR5229033,SRR5229035,SRR5229037,SRR5229039,SRR5229041,SRR5229043,SRR5229045,SRR5229047,SRR5229049,SRR5229051,SRR5229053,SRR5229055,SRR5229057,SRR5507553,SRR5507554,SRR5507555 |  sed "s/,/\n/g" | xargs -n1 -P 22 -I {} sh -c "fastq-dump {}" 


Each downloaded file needs to be demultiplexed. To do so the barcodes per sample (https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE94489&format=file&file=GSE94489%5FREADME%2Etxt) and the SRR to sample mapping needs to be provided:


.. code-block:: bash

    $ scHicDemultiplex -f "FASTQ_FILE" --srrToSampleFile samples.txt --barcodeFile GSE94489_README.txt --threads 20



After demultiplexing, each forward and reverse strand file needs to be mapped as usual in Hi-C as single-paired files:

.. code-block:: bash

    $ 

As a last step, the matrices for each cell need to be created, we use the tool 'hicBuildMatrix' from HiCExplorer:

.. code-block:: bash

    $ hicBuildMatrix

To make this step more automated, it is recommend to use either a platform like hicexplorer.usegalaxy.eu or to use a batch script:

.. code-block:: bash
    
    $ ls -1 *.bam |  xargs -n2 -P 1 -I {} sh -c "hicBuildMatrix -s {} {} --binSize 1000000 name.cool"

    