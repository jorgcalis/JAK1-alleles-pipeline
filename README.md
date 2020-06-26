# JAK1-alleles-pipeline
Pipeline to extract JAK1 allelic expression counts at single-cell resolution from an InDrops experiment.
This pipeline is part of the study described in Gruber et al. (2020 Immunity) (link to study here). 

*Step 1: Data processing of InDrops sequencing reads was performed with the indrops.py pipeline (https://github.com/indrops/indrops). From this analysis we obtained single-cell gene expression matrices, and alignment (bam) files with UMI and cell barcode attributes.

<indrops command here>

*Step 2: BAM files were processed to obtain JAK1 allele specific transcript counts with single-cell resolution. Here we use the custom build pipeline in script GITHUB_JAK1umiCountPipeline.py. For our study, we ran it following this command: 

```./GITHUB_JAK1umiCountPipeline.py
  outputDirectory=JAK1alignmentsAnalysis/ \
  sample,HealthyControl_JAK1_1,BAMalignmentFiles/PBMC_HealthyControl_JAK1_1.bam \
  sample,HealthyControl_JAK1_2,BAMalignmentFiles/PBMC_HealthyControl_JAK1_2.bam \
  sample,HealthyControl_JAK1_3,BAMalignmentFiles/PBMC_HealthyControl_JAK1_3.bam \
  sample,HealthyControl_JAK1_4,BAMalignmentFiles/PBMC_HealthyControl_JAK1_4.bam \
  sample,Patient_JAK1_1,BAMalignmentFiles/PBMC_Patient_JAK1_1.bam \
  sample,Patient_JAK1_2,BAMalignmentFiles/PBMC_Patient_JAK1_2.bam \
  sample,Patient_JAK1_3,BAMalignmentFiles/PBMC_Patient_JAK1_3.bam \
  sample,Patient_JAK1_4,BAMalignmentFiles/PBMC_Patient_JAK1_4.bam \
  > GITHUB_JAK1umiCountPipeline.out```

The input for this run are the sample names and per sample bam alignment files.
The script returns the following files for each sample:

* A per cell, per JAK1 allele count of assigned UMIs. This is the JAK1 alleles expression matrix that is used in downstream analyses. 
* An analysis on the number of identified transcripts relative to sequencing depth, to estimate if additional sequencing runs can help to increase the number of identified JAK1 transcripts.
* An overview of non-JAK1 primed genes, to facilitate the trouble-shooting of possible off-target priming.
* Statistics on the number of reads and UMIs per JAK1 allele, and on the difference between UMI sequences. For quality control purposes, and to generate run summaries for each sample. 

Results can be found in ./JAK1alignmentsAnalysis/. 
