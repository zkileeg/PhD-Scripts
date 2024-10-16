#!/bin/bash

#SBATCH --account=def-gmott
#SBATCH --cpus-per-task=1
#SBATCH --time=30:00
#SBATCH --mem=124G
#SBATCH --cpus-per-task=32
#SBATCH --array=1-2

module load seqkit

cds_file=${1}
input_seqnames_dir=${2}
output=${3}

find $input_seqname_dir -type f -name "*.sequences.txt" | sort > OGfiles.txt

cdsog=$( awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' OGfiles.txt )
#cdsog_name="${fasta_reference##*/}"

seqkit grep --threads 32 -f ${cdsog} ${cds_file} > ${cdsog##*/}.cds
