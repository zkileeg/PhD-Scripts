#!/bin/bash

#SBATCH --time=0-24:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-2

module load mmseqs2/15-6f452

gemoma_jar=${1}
fasta_reference_folder=${2}
reference_annotation=${3}
reference_genome=${4}
outdir=${5}

find $fasta_reference_folder -type f -name "*.fasta" | sort > genome_files.txt

target_genome=$( awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' genome_files.txt )
genome_name="${target_genome##*/}"



java -jar $gemoma_jar CLI GeMoMaPipeline GeMoMa.Score=ReAlign threads=32 \
	outdir=$outdir$genome_name AnnotationFinalizer.r=NO o=true \
	t=$target_genome a=$reference_annotation g=$reference_genome
