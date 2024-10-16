#!/bin/bash
#SBATCH --time=0-06:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL


module load StdEnv/2020
module load bowtie2/2.4.4

fasta_in=${1}
index_name=${2}
out_file=${3}
fastq_1=${4}
fastq_2=${5}

bowtie2-build -f $fasta_in \
/home/zkileeg/scratch/Ecotype_Data/Bowtie2_Indexes/$index_name

#Loop over the different arabidopsis ecotype files found within the accession file list.
#cat ./pair_end_sequenced/Accession_names_truncated_list.txt | while read line; do

#line="kbsmac74_illumina"


	#begin alignment with very sensitive parameters using up to 32 threads.
bowtie2 --very-sensitive -p 32 -x /home/zkileeg/scratch/Ecotype_Data/Bowtie2_Indexes/$index_name \
-S $out_file -1 $fastq_1 -2 $fastq_2


#done
