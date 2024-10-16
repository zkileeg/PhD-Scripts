#!/bin/bash

#SBATCH --time=0-24:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-26%10

module load StdEnv/2023
module load hisat2/2.2.1


fasta_reference_folder=${1}
rna_reads_dir=${2}
mapped_output_dir=${3}
fasta_index_dir=${4}

find $fasta_reference_folder -type f -name "*.fasta" | sort > genome_files.txt

fasta_reference=$( awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' genome_files.txt )
fasta_index_name="${fasta_reference##*/}"

cd $fasta_index_dir
hisat2-build -p 32 $fasta_reference $fasta_index_name

#Loop over the different arabidopsis ecotype files found within the accession file list.
#for readsfolder in $rna_reads_dir*; do

if [ ! -d $mapped_output_dir ]; then
	mkdir $mapped_output_dir
fi


	#change directory into directory containing rna reads
	cd $rna_reads_dir
	#echo "Changing into directory " $readsfolder

#get comma separated file of all fastq.gz files for RNA read files
RNAfastq_reads=$( find *.fastq | paste -s -d, - )

output_name=${mapped_output##*/}
	#begin alignment with very sensitive parameters using up to 32 threads.
hisat2 -p 32 -x $fasta_index_dir$fasta_index_name -U $RNAfastq_reads \
-S $mapped_output_dir"${fasta_index_name}_Hisat2.sam"


#done
samtools view -@ 32 -bS -F 4 $mapped_output_dir"${fasta_index_name}_Hisat2.sam" | \
samtools sort -@ 32 -o "${mapped_output_dir}${fasta_index_name}_Hisat2.bam" -
#index bam file
samtools index -@ 32 "${mapped_output_dir}${fasta_index_name}_Hisat2.bam"

rm $mapped_output_dir"${fasta_index_name}_Hisat2.sam"
