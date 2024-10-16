#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL

module load StdEnv/2023
module load apptainer/1.2.4

WORKDIR=${1}
genome_file=${2}
gene_pred=${3}
transcripts=${4}



apptainer run -C -W${SLURM_TMPDIR} -B${WORKDIR}:/home -B /scratch -B /localscratch /home/zkileeg/pasapipeline.v2.5.3.simg \
	/usr/local/src/PASApipeline/scripts/Load_Current_Gene_Annotations.dbi \
     -c alignAssembly.config -g $genome_file \
     -P $gene_pred --CPU 1

echo "Running load annotation"
#echo $working_dir
apptainer run -C -W${SLURM_TMPDIR} -B${WORKDIR}:/home -B /scratch -B /localscratch /home/zkileeg/pasapipeline.v2.5.3.simg \
	/usr/local/src/PASApipeline/Lauch_PASA_pipeline.pl \
	-c annotCompare.config -A \
	-g $genome_file \
	-t $transcripts --CPU 1
echo "running pasa pipeline"
