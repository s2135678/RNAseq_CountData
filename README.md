# RNAseq Count Data Generating Tool  

This Bash program performs quality control (FastQC) and alignment (Bowtie2) and generates a counts data from the raw RNA-sequencing data.
This tool was designed and built as part of the MSc Bioinformatics course (Bioinformatics Programming module) at University of Edinburgh (2020). 
The pipeline is implemented in 'assignment.sh'. Draft versions of the program can be found in the repository named 'draft'. 

## Overview

The program compares the RNA transcript of the parasite, _Trypanosoma brucei_, which can exist in infective and dormant forms. However, the code can be easily modified to analyse other RNA-seq data as long as the input data is in FASTQ format.    
The pipeline is mainly divided into 3 parts:

1. Quality assessment and control (FastQC)
2. Alignemnt (Bowtie2)
3. Generation of counts data 

The program outputs a text tab-delimited file containing genes on rows. 

