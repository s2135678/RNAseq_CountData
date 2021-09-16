#!/bin/bash

#Creating directories that will contain all results obtined from this analysis 
mkdir ~/Tbb
mkdir ~/Tbb/fastqc.result

#Running fastqc, results are unzipped and saved in ~/Tbb/fastqc.result  
echo "*****Currently running fastqc on all raw sequence files, please wait..." 
fastqc --extract -q -o ~/Tbb/fastqc.result /localdisk/data/BPSM/Assignment1/fastq/*.fq.gz 

#Combines summary.txt of all raw sequences into one list called summary.list
#summary.list will be the input for the while loop below  
cat  ~/Tbb/fastqc.result/*fastqc/summary.txt  > ~/Tbb/fastqc.result/summary.list


#Prints modules that raised a warning or failure from the quality assessment onto the user's screen
echo "*****Quick summary of the quality assessment. Modules that raised a warning or failure are shown below. For more information on the quality, please look at ~/Tbb/fastqc.result."
IFS=$'\t'
while read score feature filename #creating three variables 
do
if test $score == "FAIL" #if first column of summary.list contains "FAIL" prints the feature and filename to the screen 
  then 
    echo -e "${filename}: \t${feature} module has issued a failure." 
  elif test $score == "WARN" #if first column of summary.list contains "WARN" prints the feature and filename to the screen 
  then
    echo -e "${filename}: \t${feature} module has issued a warning."
  else
    continue 
fi
done < ~/Tbb/fastqc.result/summary.list  

#################################################
#################################################

#Asks user if he/she wants to continue to alignment after assessing the quality  
echo -n "Would you like to continue (y/n)?" #asks user for an input to the question 
read answer 
if [ "$answer" != "${answer#[y]}" ] #if answer is 'y' program continues  
then
  echo "*****Moving onto alignment now."
else       #if answer is 'n' program exits and returns to the terminal 
  exit 0
fi





