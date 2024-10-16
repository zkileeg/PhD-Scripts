#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL

#load modules
module load StdEnv/2020
module load bowtie2/2.4.4
module load java/17.0.2
module load pilon/1.24
module load samtools/1.16.1


#get variables from standard in arguments
genome=${1}
index_name=${2}
fastq_1=${3}
fastq_2=${4}
polished_out=${5}


#set flag to 1 so we can polish 3 times
flag=1
#Loop over the different arabidopsis ecotype files found within the accession file list.
while [ $flag -le 3 ]; do

	#set local variables to be used in the loop
        sam_file="${index_name}${flag}.sam"
	bam_file="${index_name}${flag}.bam"


	#build bowtie index. index name is one supplied by user
	bowtie2-build -f $genome "${index_name}${flag}"

	#align reads to genomes for polishing in next step
	bowtie2 --very-sensitive -p 32 -x "${index_name}${flag}" \
	-S $sam_file -1 $fastq_1 -2 $fastq_2




	#convert sam file to sorted bam file
        samtools view -@ 32 -bS -F 4 $sam_file | samtools sort -@ 32 -o $bam_file -
        #index bam file
        samtools index -@ 32 $bam_file

	#polish the genome using pilon. Give it 124G of memory. It gets angry when we use less. 

	java -Xmx124G -jar $EBROOTPILON/pilon.jar --threads 32 --genome $genome --bam $bam_file --output "${polished_out}${flag}"

	#rename for next iteration
	genome="${polished_out}${flag}.fasta"

	flag=$(( $flag + 1 ))   #iterate

done








