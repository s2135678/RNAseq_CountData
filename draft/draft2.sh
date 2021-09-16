#/bin/bash
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

#obvs al doesnt need this step so remove later
#rm -r ~/workdraft/sequence_data
#Moving reference genome and sequences, unzipped into a directory in homespace
echo "*****Moving reference genome and raw sequences into 'sequence_data' directory in your home directory. All files are unzipped."
#mkdir ~/workdraft/sequence_data
#cp /localdisk/data/BPSM/Assignment1/Tbb_genome/Tb927_genome.fasta.gz /localdisk/data/BPSM/Assignment1/fastq/*.fq.gz ~/workdraft/sequence_data
#gunzip ~/workdraft/sequence_data/Tb927_genome.fasta.gz ~/workdraft/sequence_data/*.fq.gz


#Building bowtie index, -f specifies input file as FASTA format, -q quiet 
echo "*****Currently indexing the reference genome. This may take some time..."
#bowtie2-build -q -f ~/workdraft/sequence_data/Tb927_genome.fasta Tbb_index
#mv *.bt2 ~/workdraft/sequence_data
echo "Done! The indexes can be found in 'sequence_data' directory."

#Running alignment using Bowtie2
echo "*****Currently running alignment using bowie2, this may also take some time..."
unset IFS
for i in $(ls ~/workdraft/sequence_data/*1.fq); #For all the sequence files in sequence_data directory, run:
  do
    seqname=${i::-5} #'seqname' variable defines the path to the raw sequence without last 5 chracters '_1.fq' 
    echo -e "*********\nAligning sequence ${seqname:49} to the genome..." 
    bowtie2 -x ~/workdraft/sequence_data/Tbb_index -1 ${seqname}_1.fq -2 ${seqname}_2.fq -S ${seqname:49}.sam
    samtools view -b ${seqname:49}.sam > ${seqname:49}.bam  #Converting sam files to bam, naming each one by the sequence file name
done
 
#Moving all sam and bam files into 'bowtie.result' directory in home directory 
mkdir ~/workdraft/bowtie.result
mv *.sam *.bam ~/workdraft/bowtie.result
echo "Results of alignment stored in 'bowtie_result' directory in home directory. Both in .sam and .bam format"





