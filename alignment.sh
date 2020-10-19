#/bin/bash


#Copying reference genome and sequences, moving into ~/Tbb/sequence_data and unzipping all
#echo "*****Moving reference genome and raw sequences into ~/Tbb/sequence_data. All files are unzipped."
#mkdir ~/Tbb/sequence_data
#cp /localdisk/data/BPSM/Assignment1/Tbb_genome/Tb927_genome.fasta.gz /localdisk/data/BPSM/Assignment1/fastq/*.fq.gz ~/Tbb/sequence_data
#gunzip ~/Tbb/sequence_data/Tb927_genome.fasta.gz ~/Tbb/sequence_data/*.fq.gz


#Building bowtie index, -f specifies input file as FASTA format, -q quiet 
#echo "*****Currently indexing the reference genome. This may take some time..."
#bowtie2-build -q -f ~/Tbb/sequence_data/Tb927_genome.fasta Tbb_index
#mv *.bt2 ~/Tbb/sequence_data  #moving all indexes into sequence_data directory
#echo "Done! The indexes can be found in ~/Tbb/sequence_data."

#Running alignment using Bowtie2 and convering all .sam files into .bam using samtools
#For loop runs through all sequence files and aligns to the bowtie indexes 
rm  ~/Tbb/sequence_data/*.fq.gz
echo "*****Currently running alignment using bowie2, this may also take some time..."
unset IFS
for i in $(ls ~/Tbb/sequence_data/*1.fq); #For all the sequence files in sequence_data directory, run:
  do
    seqname=${i::-5} #'seqname' variable defines the path to the raw sequence without last 5 chracters '_1.fq' 
    echo -e "*****\nAligning sequence ${seqname:43} to the genome..." 
    bowtie2 -x ~/Tbb/sequence_data/Tbb_index -1 ${seqname}_1.fq -2 ${seqname}_2.fq -S ${seqname:43}.sam
    samtools view -b ${seqname:43}.sam > ${seqname:43}.bam  #Converting sam files to bam, naming each one by the sequence file name
done
 
#Moving all .sam and .bam files into ~/Tbb/bowtie.result
mkdir ~/Tbb/bowtie.result
mv *.sam *.bam ~/Tbb/bowtie.result
echo "Results of alignment stored in ~/Tbb/bowtie.result. Both in .sam and .bam format"






