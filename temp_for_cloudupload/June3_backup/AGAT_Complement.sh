#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL


#load modules
module load StdEnv/2023
module load apptainer/1.2.4

ref_annotation=${1}
complement_annotation=${2}
output=${3}
WORKDIR=${4}

apptainer run -C -W${SLURM_TMPDIR} -B${WORKDIR}:/home -B /scratch -B /localscratch /scratch/zkileeg/Ecotype_Data/AGAT/agat_1.2.0.sif \
	agat_sp_complement_annotations.pl --ref ${ref_annotation} --add $complement_annotation -o ${output}  #complement annotations

