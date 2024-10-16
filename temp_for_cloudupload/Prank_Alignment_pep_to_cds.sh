#!/bin/bash

#SBATCH --time=0-4:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=125G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-2

module load StdEnv/2020  gcc/9.3.0
module load prank/170427

input_seq_dir=${1}
output_dir=${2}
alignment_dir=${3}

#echo $input_seq_dir
#echo $output_dir

find ${input_seq_dir} -maxdepth 1 -type f -name "*.txt.cds" | sort > cds_files.txt

input_cds=$( awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' cds_files.txt )


find ${alignment_dir} -maxdepth 1 -type f -name "*.txt.cds" | sort > pep_aligned_files.txt

aligned_in=$( awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' pep_aligned_files.txt )


output_seq="${output_dir}${input_cds##*/}_aligned.codon"

echo $input_cds
echo $output_seq

#prank -d=$input_cds -o=$output_seq -codon -F
prank -convert -d=$alignment_in -dna=$input_cds -o=$output_seq
