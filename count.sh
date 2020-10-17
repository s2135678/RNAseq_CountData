#!/bin/bash

#sorting sam and bam files into slender and stumpy 
rm -r ~/workdraft/bowtie.result
rm -r ~/workdraft/final.result
mkdir {~/workdraft/bowtie.result/slender,~/workdraft/bowtie.result/stumpy}
for z in $(ls ~/workdraft/bowtie.result/*am);
do
  num=${z:49:3}
  if test $num -ge "220"
  then
    mv $z ~/workdraft/bowtie.result/stumpy
  else
    mv $z ~/workdraft/bowtie.result/slender
  fi
done

#Generating counts data using bedtools 
mkdir ~/workdraft/final.result
bedtools intersect -a /localdisk/data/BPSM/Assignment1/Tbbgenes.bed -b ~/workdraft/bowtie.result/slender/*.bam -c > ~/workdraft/final.result/slender.bam
bedtools intersect -a /localdisk/data/BPSM/Assignment1/Tbbgenes.bed -b ~/workdraft/bowtie.result/stumpy/*.bam -c > ~/workdraft/final.result/stumpy.bam

