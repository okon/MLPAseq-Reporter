MLPA-seq Reporter
============

## 1. About
MLPA-seq reporter is a tool designed for analysis of MLPA-seq data.

**Author:**  
Olga Kondrashova  
*Department of Pathology, The University of Melbourne.*


## 2. Requirements
+ **System Requirements:**
* Linux operating system with BASH shell (for AmpliVar)

+ **Software Requirements:**
* AmpliVar
* R version 3.1.2
* R packages:
- ggplot2
- reshape
- getopt
- matrixStats

## 3. Usage and Program Options
+ **File naming and formats**
* Paired R1 and R2 FASTQ files in the input directory must have the following format:
RunDate_InstrumentN_RunN_LibraryID_SampleID-ProbeMix_ExperimentID_Assay_LaneN_ReadN.fastq.gz
For example:
150101_M00001_0001_LB001_S01-MixA_E01_MLPAseq_L001_R1.fastq.gz
150101_M00001_0001_LB001_S01-MixA_E01_MLPAseq_L001_R2.fastq.gz
* Genotype file is a four-column, tab-separated text file with the first three columns describing the 
probe pair and the fourth column is the corresponding sequence being genotyped (matched). A generic genotype file example is provided in AmpliVar instructions. A typical genotype file required for MLPA-seq Reporter analysis is described here. The probe-pair description provided in the first column is critical for proper script performance - it lists the probe number (arbitrary), gene name, exon/intron number (if available), the probe number for that region (a,b,c... or 1,2,3...), adapter used (Nilmn - Nextera Illumina), probe mix (MixA - detection probe mix or MixB - confirmation probe mix), and whether the probe pair is control (ctrl) or not (nothing). Each value is separated by an underscore. Second column is an arbitrary letter (usually M for middle of the probes). Third column is the probe-pair position. Forth column is the reference matching sequence (generally the proximal 8 bases of left and right probe). For example (without the header row):

	   | Column 1         | Column 2     | Column 3   | Matching sequence               |   
	   | ---------------- | ------------ | ---------- | ------------------------------- |  
       | 001_CFTR_a_Nilmn_MixAMixB_ctrl  | M          | chr7:117307078-117307150       | TGCTCTGAAAGAGGAG |  
       | 014_BRCA1_ex2_a_Nilmn_MixA_ctrl | M          | chr17:41276060-41276123      | AACGCGAAGAGCAGAT |  
       | 015_BRCA1_ex2_b_Nilmn_MixB_ctrl | M          | chr17:41276019-41276094      | AAATCTTAGAGTGTCC |  
       | 151_CCNE1_ex11_1_Nilmn_MixAMixB | M          | chr19:30313375-30313438      | CCATGGCAAATGGAAC |  
       | 152_CCNE1_ex6_1_Nilmn_MixAMixB  | M          | chr19:30308310-30308385      | ACATGATTTTCCAGAC |   
* Flanks file
