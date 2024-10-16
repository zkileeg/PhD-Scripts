#!/bin/bash
#SBATCH --time=0-24:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-44

module load StdEnv/2020
#module load stringtie/2.2.1

bam_dir_in=${1}
gtf_output_dir=${2}

#find all instances of files with .bam extension and output it to a list
if [[ $SLURM_ARRAY_TASK_ID -eq 1 ]]; then
	find ${bam_dir_in} -maxdepth 1 -type f -name "*.bam" | sort > bam_files.txt
fi
#using array job and the internal slurm id counter ($SLURM_ARRAY_TASK_ID = the current array task number),
#set the bam_reference path to the line found at that point in the bam_files text file
bam_reference=$( awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' bam_files.txt )
#extract only the name and extension of the bam file so we a hve a name for output
bam_out=${gtf_output_dir}"${bam_reference##*/}.gtf"

echo $bam_reference
echo $bam_out

#assemble transcripts using stringtie, with the input as the bam reference. Assemble the transcripts de novo
=/home/zkileeg/stringtie-2.2.2.Linux_x86_64/stringtie -o ${bam_out} -p 32 $bam_reference
#stringtie -o ${bam_out} -p 32 $bam_reference
