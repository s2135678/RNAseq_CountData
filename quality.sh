#!/bin/bash

#Remvoing old and creating folder that will contain the fastqc result
#No need to remove for Al, he wouldn't have it to begin with 
rm -r fastqc.result
mkdir fastqc.result

#Running fastqc, results are saved in the 'fastqc.result' directory in home directory
echo "*****Currently running fastqc on all raw sequence files, please wait..." 
fastqc --extract -q -o ~/workdraft/fastqc.result /localdisk/data/BPSM/Assignment1/fastq/*.fq.gz 

#Combines summary.txt of all raw sequences into one list
cat  ~/workdraft/fastqc.result/*fastqc/summary.txt  > ~/workdraft/fastqc.result/summary.list

#########################

echo "*****Quick summary of quality assessment. Analysis modules that raised a warning or failure are shown below. For more information on the quality, please look at 'fastqc.result' directory in your home directory"

#Prints analysis mododules that raised a warning or failure onto the screen
IFS=$'\t'
while read score feature filename
do
if test $score == "FAIL"
  then 
    echo -e "${filename}: \t${feature} module has issued a failure"
  elif test $score == "WARN"
  then
    echo -e "${filename}: \t${feature} module has issued a warning"
  else
    continue 
fi
done < ~/workdraft/fastqc.result/summary.list  ##change here to home directory

#Asks user if he/she wants to continue to alignment after assessing the quality 
#If answer is 'y' program continues 
#If answer is 'n' it exits the shell script and returns to the terminal 
echo -n "Would you like to continue (y/n)?"
read answer
if [ "$answer" != "${answer#[Yy]}" ] 
then
  echo "Moving onto the alignment now"
else
  exit 0
fi





