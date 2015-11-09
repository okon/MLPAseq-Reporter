MLPA-seq Reporter
=================


## 1. About
MLPA-seq reporter is a tool designed for analysis of MLPA-seq data.

**Author:**
Olga Kondrashova
*Department of Pathology, The University of Melbourne.*


## 2. Requirements
+ **System Requirements:**
 * Linux operating system with BASH shell (for AmpliVar)
 * R version 3.1.2

+ **Software Dependencies:**
 * AmpliVar


## 3. Usage and Program Options
+ **File naming and formats**
	* Paired R1 and R2 FASTQ files in the input directory must have the following format:
		RunDate_InstrumentN_RunN_LibraryID_SampleID-ProbeMix_ExperimentID_Assay_LaneN_ReadN.fastq.gz
	  For example:
		150101_M00001_0001_LB001_S01-MixA_E01_MLPAseq_L001_R1.fastq.gz
		150101_M00001_0001_LB001_S01-MixA_E01_MLPAseq_L001_R2.fastq.gz
	* Genotype file is a four-column, tab-separated text file with the first three columns describing the 
	   probe pair and the fourth column is the corresponding sequence being genotyped (matched). The probe-pair description provided in the first column is critical for proper normalisation, it lists the probe number (arbitrary), gene name, exon number, the probe number in that exon (usually a,b,c etc.), probe adapter (Nilmn - Nextera Illumina), probe mix (MixA - detection probe mix or MixB - confirmation probe mix), and whether the probe pair is control (ctrl) or not (nothing). Each value is separated by underscore. Second column is arbitrary letter (usually M for middle of the probes). Third column is the probe pair position. Forth column
generally the proximal 10 bp of each probe. For example (without header row):
