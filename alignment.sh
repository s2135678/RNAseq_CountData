#/bin/bash

#obvs al doesnt need this step so remove later
#rm -r ~/workdraft/sequence_data
#Moving reference genome and sequences, unzipped into a directory in homespace
echo "*****Moving reference genome and sequences into sequence_data directory on your home directory. All files are unzipped"
#mkdir ~/workdraft/sequence_data
#cp /localdisk/data/BPSM/Assignment1/Tbb_genome/Tb927_genome.fasta.gz /localdisk/data/BPSM/Assignment1/fastq/*.fq.gz ~/workdraft/sequence_data
#gunzip ~/workdraft/sequence_data/Tb927_genome.fasta.gz ~/workdraft/sequence_data/*.fq.gz
echo "Done!"

##### In the real script you have to make index 
#Building bowtie index, -f specifies input file as FASTA format, -q quiet 
echo "*****Currently indexing the reference genome. This may take some time..."
#bowtie2-build -q -f ~/workdraft/sequence_data/Tb927_genome.fasta Tbb_index
#mv *.bt2 ~/workdraft/sequence_data
echo "Done! The indexes can be found in sequence_data directory."

#Running alignment using Bowtie2
echo "*****Currently running alignment using bowie2, this also may take some time..."
bowtie2 -x ~/workdraft/sequence_data/Tbb_index -1 ~/workdraft/sequence_data/216_L8_1.fq -2 ~/workdraft/sequence_data/216_L8_2.fq -S eg1.sam

#Converting sam into bam
echo "*****Converting sam files into bam."
mkdir ~/workdraft/bowtie.result
samtools view -b eg1.sam > eg1.bam
mv *.sam *.bam ~/workdraft/bowtie.result
echo "Results of alignment stored in bowtie_result directory in home directory both in .sam and .bam format"