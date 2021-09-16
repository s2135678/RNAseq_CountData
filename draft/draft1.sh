#!/bin/bash

###Remvoing old and creating folder that will contain the fastqc result
#No need to remove for Al, he wouldn't have it to begin with 
rm -r fastqc.result
mkdir fastqc.result

###Running fastqc, --extract to unzip the files and -q for no message output
echo "*****Currently running fastqc on all raw sequence files, please wait a few moments*****" ##change here to create fastqc.result in home directory####
fastqc --extract -q -o ~/workdraft/fastqc.result /localdisk/data/BPSM/Assignment1/fastq/*.fq.gz 

###Moves basic summary.txt files of all rae sequences into one text file 
###The text file is moved into fastqc.result directory to prevent overcrowding the working directory 
##change here to access fastqc.result file in home directory## 
cat  ~/workdraft/fastqc.result/*fastqc/summary.txt  > summary.list
mv summary.list ~/workdraft/fastqc.result

#########################
#########################

echo "*****Quick summary of the quality test for raw sequences*****"
echo "*****For more detailed information of each sequence file, look at fastqc.result directory in your home directory*****"

#Print to the screen if some features scored low in quality assessment 
IFS=$'\t'
while read score feature filename
do
if test $score == "FAIL"
  then 
    echo -e "${filename} has failed the quality assessment  ${feature}"
  elif test $score == "WARN"
  then
    echo -e "${filename} has been warned for the quality assessment  ${feature}"
  else
    continue 
fi
done < ~/workdraft/fastqc.result/summary.list  ##change here to home directory

#Asks user if he/she wants to continue to alignment after assessing the quality 
#read -p "Would you like to continue (y/n)?" answer
#if $answer == "y"
#  then 
#   echo "Moving onto the alignment now"
   
#   else
#    kill
#fi
  
####################################
####################################
#Alignment   
  
#obvs al doesnt need this step so remove later
rm -r ~/workdraft/sequence_data
#Moving reference genome and sequences, unzipped into a directory in homespace
echo "Moving reference genome and sequences into sequence_data directory on your home directory. All files are unzipped"
mkdir ~/workdraft/sequence_data
cp /localdisk/data/BPSM/Assignment1/Tbb_genome/Tb927_genome.fasta.gz /localdisk/data/BPSM/Assignment1/fastq/*.fq.gz ~/workdraft/sequence_data
gunzip ~/workdraft/sequence_data/Tb927_genome.fasta.gz ~/workdraft/sequence_data/*.fq.gz
echo "Done!"

##### In the real script you have to make index 
#Building bowtie index, -f specifies input file as FASTA format, -q quiet 
echo "Currently indexing the reference genome. This may take some time..."
bowtie2-build -q -f ~/workdraft/sequence_data/Tb927_genome.fasta Tbb_index
mv *.bt2 ~/workdraft/sequence_data
echo "Done! The indexes can be found in sequence_data directory."

#Running alignment using Bowtie2
echo "Currently running alignment using bowie2, this also may take some time..."
bowtie2 -x ~/workdraft/sequence_data/Tbb_index -1 ~/workdraft/sequence_data/216_L8_1.fq -2 ~/workdraft/sequence_data/216_L8_2.fq -S eg1.sam

#Converting sam into bam
samtools view -b eg1.sam > eg1.bam
mv *.sam *.bam ~/workdraft/bowtie.result
echo "Results of alignment stored in bowtie_result directory in home directory both in .sam and .bam format"    






