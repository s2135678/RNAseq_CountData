#!/bin/bash

###Remvoing old and creating folder that will contain the fastqc result
rm -r fastqc.result
mkdir fastqc.result

###Running fastqc, --extract to unzip the files and -q for no message output
echo "*****Currently running fastqc on all raw sequence files, please wait a few moments*****" 
fastqc --extract -q -o ~/workdraft/fastqc.result /localdisk/data/BPSM/Assignment1/fastq/*.fq.gz 

###Moves basic summary.txt files of all rae sequences into one text file 
###The text file is moved into fastqc.result directory to prevent overcrowding the working directory 
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
done < ~/workdraft/fastqc.result/summary.list

#Asks user if he/she wants to continue to alignment after assessing the quality 
read -p "Would you like to continue (y/n)?" answer
if $answer == "y"
  then 
   echo "Moving onto the alignment now"
   else
    kill
fi
  
    






