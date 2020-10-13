#!/bin/bash

rm -r fastqc.result
mkdir fastqc.result
fastqc --extract -o ~/Assignment1_draft/fastqc.result /localdisk/data/BPSM/Assignment1/fastq/216_L8_1.fq.gz 

cat  ~/Assignment1_draft/fastqc.result/216_L8_1_fastqc/summary.txt  > summary.list
echo "Pointing out the features that fail or got warn"

IFS=$'\t'
while read score feature filename
do
if test $score == "FAIL"
  then 
    echo -e "${filename} has fail the quality assessment  ${feature}"
  elif test $score == "WARN"
  then
    echo -e "${filename} has warn  the quality assessment  ${feature}"
  else
    continue 
fi
done < summary.list

