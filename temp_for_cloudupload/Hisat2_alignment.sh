#!/bin/bash

#SBATCH --time=0-24:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=125G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020
module load hisat2/2.2.1


fasta_reference=${1}
fasta_index_name=${2}
rna_reads_dir=${3}
mapped_output=${4}
fasta_index_dir=${5}

cd $fasta_index_dir
hisat2-build -p 32 $fasta_reference $fasta_index_name

#Loop over the different arabidopsis ecotype files found within the accession file list.
#for readsfolder in $rna_reads_dir*; do

if [ ! -d $mapped_output ]; then
	mkdir $mapped_output
fi


	#change directory into directory containing rna reads
	cd $rna_reads_dir
	#echo "Changing into directory " $readsfolder

#get comma separated file of all fastq.gz files for RNA read files
RNAfastq_reads=$( find *.fastq | paste -s -d, - )

output_name=${mapped_output##*/}
	#begin alignment with very sensitive parameters using up to 32 threads.
hisat2 -p 32 -x $fasta_index_dir$fasta_index_name -U $RNAfastq_reads \
-S $mapped_output/"${fasta_index_name}_Hisat2.sam"


#done
