#!/bin/bash

###########################################
###### Quality assessment #################
###########################################

#Creating directories that will contain all results obtined from this analysis in the current workding directory
mkdir ./Tbb
mkdir ./Tbb/fastqc.result

#Running fastqc, results are unzipped and saved in ./Tbb/fastqc.result  
echo "*****Currently running fastqc on all raw sequence files, please wait..." 
fastqc --extract -q -o ./Tbb/fastqc.result /localdisk/data/BPSM/Assignment1/fastq/*.fq.gz 

#Combines summary.txt of all raw sequences into one list called summary.list
#summary.list will be the input for the while loop below  
cat  ./Tbb/fastqc.result/*fastqc/summary.txt  > ./Tbb/fastqc.result/summary.list


#Prints modules that raised a warning or failure from the quality assessment onto the user's screen
echo "*****Quick summary of the quality assessment. Modules that raised a warning or failure are shown below. For more information on the quality, please look at ./Tbb/fastqc.result."
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
done < ./Tbb/fastqc.result/summary.list  

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

#################################################
########### Alignment ###########################
#################################################


#Copying reference genome and sequences, moving into ./Tbb/sequence_data and unzipping all
echo "*****Moving reference genome and raw sequences into ./Tbb/sequence_data. All files are unzipped."
mkdir ./Tbb/sequence_data
cp /localdisk/data/BPSM/Assignment1/Tbb_genome/Tb927_genome.fasta.gz /localdisk/data/BPSM/Assignment1/fastq/*.fq.gz ./Tbb/sequence_data
gunzip ./Tbb/sequence_data/Tb927_genome.fasta.gz ./Tbb/sequence_data/*.fq.gz


#Building bowtie index, -f specifies input file as FASTA format, -q quiet 
echo "*****Currently indexing the reference genome. This may take some time..."
bowtie2-build -q -f ./Tbb/sequence_data/Tb927_genome.fasta Tbb_index
mv *.bt2 ./Tbb/sequence_data  #moving all indexes into sequence_data directory
echo "Done! The indexes can be found in ./Tbb/sequence_data."

#Running alignment using Bowtie2 and convering all .sam files into .bam using samtools
#For loop runs through all sequence files and aligns to the bowtie indexes 
echo "*****Currently running alignment using bowtie2, this may also take some time..."
unset IFS
for i in $(ls ./Tbb/sequence_data/*1.fq); #For all the sequence files in sequence_data directory, run:
  do
    seqname=${i::-5} #'seqname' variable defines the path to the raw sequence without the last 5 chracters '_1.fq' 
    echo -e "*****\nAligning sequence ${seqname:20} to the genome..." 
    bowtie2 -x ./Tbb/sequence_data/Tbb_index -1 ${seqname}_1.fq -2 ${seqname}_2.fq -S ${seqname:20}.sam
    samtools view -b ${seqname:20}.sam > ${seqname:20}.bam  #Converting sam files to bam, naming each one by the sequence file name
done
 
#Moving all .sam and .bam files into ./Tbb/bowtie.result
mkdir ./Tbb/bowtie.result
mv *.sam *.bam ./Tbb/bowtie.result
echo "Results of alignment stored in ./Tbb/bowtie.result. Both in .sam and .bam format"

##################################################
########## Counting ##############################
##################################################

#Sorting the .sam and .bam files into 'slender' and 'stumpy' directories within bowtie.result
mkdir {./Tbb/bowtie.result/slender,./Tbb/bowtie.result/stumpy} #makes 'slender' 'stumpy' directories
for z in $(ls ./Tbb/bowtie.result/*am); #z variable defines path to all .sam and .bam files in ./Tbb/bowtie.result
do
  num=${z:20:3}  #num variable contains just the three numbers of the sequence file names eg. 216
  if test $num -ge "220"   
  then
    mv $z ./Tbb/bowtie.result/stumpy  #if num is greater or equal to 220, the files are moved into 'stumpy' directory
  else
    mv $z ./Tbb/bowtie.result/slender #any other sequence files are moved into 'slender' directory 
  fi
done

#running bedtools to count number of reads that align to the genome that code for genes 
#bedtools intersect is run separately on all .bam files in 'slender' and 'stumpy' directories  
echo "*****Moving onto the final stage: counting the number of reads that align to the genome..." 
mkdir ./Tbb/final.result #creating directory, which will contain all final results
bedtools intersect -a /localdisk/data/BPSM/Assignment1/Tbbgenes.bed -b ./Tbb/bowtie.result/stumpy/*.bam -c > ./Tbb/final.result/stumpy.bam #running for all .bam files in 'stumpy directory'
bedtools intersect -a /localdisk/data/BPSM/Assignment1/Tbbgenes.bed -b ./Tbb/bowtie.result/slender/*.bam -c > ./Tbb/final.result/slender.bam #running for all .bam files in 'slender directory'

#Producing one text files containing the gene name, average count for slender and average count for stumpy
#From stumpy.bam, each gene name (field 4) and the corresponding count (field 7) are stored in an array called stump
#If the gene name (field 4) of slender.bam matches the gene names stored in the array stump, it prints the gene name, slender gene count /3 and stumpy gene count/3
#The results are stored in ./Tbb/final.result/result.txt 
awk -F'\t' 'NR==FNR{stump[$4]=$7;next}stump[$4]{print $4"\t"$7/3"\t"stump[$4]/3}' ./Tbb/final.result/stumpy.bam ./Tbb/final.result/slender.bam > ./Tbb/final.result/result.txt  
echo -e "gene_name\tslender_mean_count\tstumpy_mean_count" > ./Tbb/final.result/FINAL.txt #adding header lines to the count data
cat ./Tbb/final.result/result.txt >> ./Tbb/final.result/FINAL.txt
rm ./Tbb/final.result/result.txt #removing the text fle without header

echo "*****Everything is completed. The text file containing the count data is called 'FINAL.txt' and it can be found in ./Tbb/final.result.Thank you!"
