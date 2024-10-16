#!/bin/bash
#SBATCH --time=0-12:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=164G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-4

module load StdEnv/2020  gcc/9.3.0
module load exonerate/2.4.0

#get input from cmd line
protein_reference=${1}
#target_assembly_folder=${2}
#output_folder=${3}
redo_file=${2}

#if [ "${output_folder: -1}" != "/" ]; then
#	output_folder="${output_folder}/"
#fi

#$find ${target_assembly_folder} -type f -name *.fasta -maxdepth 1 | sort > genome_files.txt
#SLURM_ARRAY_TASK_ID=1

target_assembly=$(awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' ${redo_file} | sed "s/_exonerate.gff//g" )
#target_assembly=$(echo ${target_assembly} | sed "s/_exonerate.gff//g")

#output_folder="${target_assembly%/*}/"
#echo $protein_reference
echo $target_assembly
echo $output_folder
echo ${target_assembly}_exonerate.gff
#output_name="${target_assembly##*/}_exonerate.gff"

exonerate -m protein2genome -c 32 --querytype protein --targettype dna --percent 70 --minintron 10 --maxintron 60000 --showtargetgff TRUE --showvulgar no --showalignment no \
 $protein_reference $target_assembly > "${target_assembly}_exonerate.gff"
