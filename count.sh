#!/bin/bash

#Sorting the .sam and .bam files into 'slender' and 'stumpy' directories within bowtie.result
mkdir {~/Tbb/bowtie.result/slender,~/Tbb/bowtie.result/stumpy} #makes 'slender' 'stumpy' directories
for z in $(ls ~/Tbb/bowtie.result/*am); #z variable defines path to all .sam and .bam files in ~/Tbb/bowtie.result
do
  num=${z:43:3}  #num variable contains just the three numbers of the sequence file names eg. 216
  if test $num -ge "220"   
  then
    mv $z ~/Tbb/bowtie.result/stumpy  #if num is greater or equal to 220, the files are moved into 'stumpy' directory
  else
    mv $z ~/Tbb/bowtie.result/slender #any other sequence files are moved into 'slender' directory 
  fi
done

#running bedtools to count number of reads that align to the genome that code for genes 
#bedtools intersect is run separately on all .bam files in 'slender' and 'stumpy' directories  
echo "*****Moving onto the final stage: counting the number of reads that align to the genome..." 
mkdir ~/Tbb/final.result #creating directory, which will contain all final results data  
bedtools intersect -a /localdisk/data/BPSM/Assignment1/Tbbgenes.bed -b ~/Tbb/bowtie.result/stumpy/*.bam -c > ~/Tbb/final.result/stumpy.bam #running for all .bam files in 'stumpy directory'
bedtools intersect -a /localdisk/data/BPSM/Assignment1/Tbbgenes.bed -b ~/Tbb/bowtie.result/slender/*.bam -c > ~/Tbb/final.result/slender.bam #running for all .bam files in 'slender directory'

#Producing one text files containing the gene name, average count for slender and average count for stumpy

awk -F'\t' 'NR==FNR{stump_gene[$4]=$7;next}stump_gene[$4]{print $4"\t"$7/3"\t"stump_gene[$4]/3}' ~/Tbb/final.result/stumpy.bam ~/Tbb/final.result/slender.bam > ~/Tbb/final.result/result.txt
echo "*****Everything is completed. The text file containing the count data is called 'result.txt' and it can be found in ~/Tbb/final.result.Thank you!"

