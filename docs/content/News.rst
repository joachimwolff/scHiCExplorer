News
====

** 5 August 2020**

Release of version 5:

- Better clustering with MinHash: More accurate, faster loading times, less memory, additional PCA
- In general: Fast loading of matrices. This version is up to 19 times faster (10 kb matrices scHiCExplorer version 4: 58 minutes. Version 5: < 3 minutes)
- Support for scool format as defined with version 0.8.9 of cooler
- scHicManageScool: Tool to update old scool (version 4 or less of scHiCExplorer) to new version. Option to extract a matrix to a single cool
- scHicConvertFormat: Tool to convert a scool matrix to the file and folder structure scHiCluster needs
- Option to plot PC1 vs PC2 on cluster results for scHicCluster and scHicClusterMinHash
- scHicClusterMinHash: New additional cluster algorithms: birch, agglomerative (ward, average, single, complete)
- scHicCluster and scHicClusterMinHash: Option to work on intra-chromosomal data only
- Multiple bug fixes
- Improved plotting of scHicPlotConsensusMatrices and scHicPlotClusterProfile

** 8 March 2020**

Release of version 4:

- Fixing a bug in scHicDemultiplex
- Improving the documentation on how to download the FASTQ files for example.

**7th March 2020**

Release of version 3:

- Change datatype name mcool to scool. 
- Change tool scHicMergeToMCool to scHicMergeToSCool
- Add tool scHicCreateBulkMatrix to create the bulk matrix out of the individual single-cell matrices

**24th February 2020**

Release of version 2.

- Option to define how many nearest neighbors in scHicCluster or scHicClusterMinHash should be computed
- scHicClusterMinHash: Partially loading of data from Python to C++ to decrease memory footprint to a 1/4. 

**30th November 2019**

Release of version 1 of the single-cell HiCExplorer.


**06th November 2019**

Creation of the single-cell HiCExplorer documentation.