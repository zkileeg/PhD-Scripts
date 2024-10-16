#!/bin/bash
#SBATCH --time=0-24:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-44


module load StdEnv/2020
module load samtools/1.16.1


bam_dir_SE=${1}
bam_dir_PE=${2}
output_dir=${3}

#Loop over the different arabidopsis ecotype files found within the accession file list.

if [[ $SLURM_ARRAY_TASK_ID -eq 1 ]]; then

	find ${bam_dir_SE} -type f -maxdepth 1 -name "*.bam"  | sort > bam_files_SE.txt
	find ${bam_dir_PE} -type f -maxdepth 1 -name "*.bam"  | sort > bam_files_PE.txt

fi	#find all instances of files with .bam extension and output it to a list




#using array job and the internal slurm id counter ($SLURM_ARRAY_TASK_ID = the current array task number),
#set the bam_reference path to the line found at that point in the bam_files text file
bam_file_SE=$( awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' bam_files_SE.txt )
bam_file_PE=$( awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' bam_files_PE.txt )
#echo ${bam_files[@]}
#convert sam file to sorted bam file
output_name="${bam_file_SE##*/}_Merged.bam"

samtools merge -@ 32 -f -o ${output_dir}${output_name} $bam_file_SE $bam_file_PE
#index bam file
samtools index -@ 32 ${output_dir}${output_name}



