MLPA-seq Reporter
============

## 1. About
MLPA-seq reporter (aka MLPAbrary_analysis) is a tool designed for analysis of MLPA-seq data. This tool performs batch analysis of a single experiment (single sequencing run), by normalising each MLPA-seq probe pair by depth within each sample and by control samples. This tool produces a summary report for each sample that lists mean ratio and standard deviation for each analysed exon (or target region), as defined in the provided genotype file, it also produces a graphical output that visualises the same data, and a raw coverage file for all samples run in that experiment/sequencing run.

The detailed description of the normalisation process is described in
[publication link](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0143006).

If you use the MLPA-seq reporter, please cite the mentioned paper.  
Kondrashova O, Love CJ, Lunke S, Hsu AL, Australian Ovarian Cancer Study (AOCS) Group, Waring PM, et al. (2015) High-Throughput Amplicon-Based Copy Number Detection of 11 Genes in Formalin-Fixed Paraffin-Embedded Ovarian Tumour Samples by MLPA-Seq. PLoS ONE 10(11): e0143006. 

**Author:**  
Olga Kondrashova :girl:  
*Department of Pathology, The University of Melbourne.*


## 2. Requirements
**System Requirements:**  

* Linux operating system with BASH shell (for AmpliVar)

**Software Requirements:**  

* AmpliVar  
* R version 3.1.2  
* R packages:  
	- ggplot2
	- getopt
	- matrixStats

## 3. Usage and Program Options
**File naming and formats**  

* Paired R1 and R2 FASTQ files in the input directory must have the following format:
RunDate_InstrumentN_RunN_LibraryID_SampleID-ProbeMix_ExperimentID_Assay_LaneN_ReadN.fastq.gz
For example:
150101_M00001_0001_LB001_S01-MixA_E01_MLPAseq_L001_R1.fastq.gz  
150101_M00001_0001_LB001_S01-MixA_E01_MLPAseq_L001_R2.fastq.gz  
150101_M00001_0001_LB001_S01-MixB_E01_MLPAseq_L001_R1.fastq.gz  
150101_M00001_0001_LB001_S01-MixB_E01_MLPAseq_L001_R2.fastq.gz  
  
	NTC samples need to have 'NTC' in place of SampleID field, while positive control samples used for ratio normalisation need to have 'control' tag in place of Sample ID field.

* Genotype file (aka usual suspects file) is a four-column, tab-separated text file with the first three columns describing the 
probe pair and the fourth column is the corresponding sequence being genotyped (matched). A generic genotype file example is provided in AmpliVar instructions. A typical genotype file required for MLPA-seq Reporter analysis is described here. The probe-pair description provided in the first column is critical for proper script performance - it lists the probe number (arbitrary), gene name, exon/intron number (if available), the probe number for that region (a,b,c... or 1,2,3...), adapter used (Nilmn - Nextera Illumina), probe mix (MixA - detection probe mix or MixB - confirmation probe mix), and whether the probe pair is control (ctrl) or not (nothing). Each value is separated by an underscore. Second column is an arbitrary letter (usually M for middle of the probes). Third column is the probe-pair position. Forth column is the reference matching sequence (generally the proximal 8 bases of left and right probe). For example (without the header row):      

    | Column 1                          |  Column 2  |  Column 3                  |  Variant sequence  | 
    |-----------------------------------|------------|----------------------------|--------------------| 
    |  001_CFTR_a_Nilmn_MixAMixB_ctrl   |  M         |  chr7:117307078-117307150  |  TGCTCTGAAAGAGGAG  | 
    |  014_BRCA1_ex2_a_Nilmn_MixA_ctrl  |  M         |  chr17:41276060-41276123   |  AACGCGAAGAGCAGAT  | 
    |  015_BRCA1_ex2_b_Nilmn_MixB_ctrl  |  M         |  chr17:41276019-41276094   |  AAATCTTAGAGTGTCC  | 
    |  151_CCNE1_ex11_1_Nilmn_MixAMixB  |  M         |  chr19:30313375-30313438   |  CCATGGCAAATGGAAC  | 
    |  152_CCNE1_ex6_1_Nilmn_MixAMixB   |  M         |  chr19:30308310-30308385   |  ACATGATTTTCCAGAC  |  

* Primer sequence file (probes/flanks file) must take the following form in tab-separated text as shown in the table below (without the header row):  

    | Amplicon                          |  Length  |  Coordinate                  |  Flanking primer sequence      | 
    |-----------------------------------|----------|------------------------------|--------------------------------| 
    |  001_CFTR_a_Nilmn_MixAMixB_ctrl   | 72       |  chr7:117307078-117307150    |  (GAACTCAAGCAA.*GAGGTGCAAGAG)  | 
    |  014_BRCA1_ex2_a_Nilmn_MixA_ctrl  | 63       |  chr17:41276060-41276123     |  (CATAGCATTAAT.*ATTTCTTTCTGT)  | 
    |  015_BRCA1_ex2_b_Nilmn_MixB_ctrl  | 75       |  chr17:41276019-41276094     |  (GCGTTGAAGAAG.*AGTCAGCACAAG)  | 
    |  151_CCNE1_ex11_1_Nilmn_MixAMixB  | 63       |  chr19:30313375-30313438     |  (CAGTTTTGAGCT.*TTTGCCCAGCTA)  | 
    |  152_CCNE1_ex6_1_Nilmn_MixAMixB   | 75       |  chr19:30308310-30308385     |  (AAAGTGCTGATC.*TTGACACAGTTC)  | 
        

**Usage**  

To install MLPA-seq Reporter, simply download the package:
         
    git clone https://github.com/okon/MLPAseq-Reporter.git

* Make sure AmpliVar is installed and is runnning properly by following [AmpliVar](https://github.com/alhsu/AmpliVar) instructions.   

* Run AmpliVar Genotyping mode on the FASTQ files required for analysis.  

        [/PATH/TO/AMPLIVAR]/bin/universal/amplivar_wrapper.sh -m GENOTYPING \
        --input [INPUTDIR] \
        --output [OUTPUTDIR] \
        --probes [PROBES] \
        --suspects  [SUSPECTS] \
        -d NEXTERA    

* Move all genotypes files (files are located in a genotype subfolder with file ending - _merged_seqprep.fastq.fna_1_num_grp_Genotypes.txt) to one directory.  

        MLPAbrary_analysis.R [-[-help|h]] [-[-InputDir|I] <character>] [-[-ControlDir|C] <character>] 
        [-[-OutputDir|O] <character>] [-[-Suspects|S] <character>] [-[-verbose|v]] [-[-Version|V]]  

**Options** 
  
    -h|--help          Print this helpful helpscreen  
    -I|--InputDir      Required: input directory of AmpliVar analysis with genotype files  
    -C|--ControlDir    Required: directory with control samples - genotype files  
    -O|--OutputDir     Required: output directory for summary reports  
    -S|--Suspects      Usual suspects file used for MLPAbrary analysis  
    -v|--verbose       Print progress messages  
    -V|--Version       Display script Version  
**Test run** 

Since FASTQ files and the amplivar directory could not be provided due to the large file size, these steps have already been performed:  
       
     mkdir ./amplivar
     
    [/PATH/TO/AMPLIVAR]/bin/universal/amplivar_wrapper.sh -m GENOTYPING \  
    --input [/PATH/TO/FASTQ] \
    --output ./amplivar \
    --probes MLPAbrary_big_flanks.txt \
    --suspects MLPAbrary_middle_long_suspects.txt \
    -d NEXTERA 
    
	mv ./amplivar/*/genotype/*Genotypes.txt ./MLPAseq_run_genotypes

Run MLPAbrary_analysis script:  

    MLPAbrary_analysis_0.6.1.R \
    -I MLPAseq_run_genotypes \
    -O Test_results \
    -S MLPAbrary_middle_long_suspects.txt
    
Check if the output files in the Test_results folder are the same as the files in the Ref_results folder:

		for T in ./Test_results/*.txt ; do   
			B=`basename $T`;  
			diff ./Ref_results/$B $R;     
		done  
 
      
       
