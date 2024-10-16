#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL

module load StdEnv/2023
module load apptainer/1.2.4

WORKDIR=${1}
sample_name=${2}
genome_file=${3}
weights_file=${4}
gene_pred=${5}
transcripts=${6}
proteins=${7}


#echo $working_dir
apptainer run -C -W${SLURM_TMPDIR} -B${WORKDIR}:/home -B /scratch -B /localscratch /home/zkileeg/EVidenceModeler.v2.1.0.simg EVidenceModeler \
                   --exec_dir $WORKDIR \
                   --sample_id $sample_name \
                   --genome $genome_file \
                   --weights $weights_file \
                   --gene_predictions $gene_pred \
                   --transcript_alignments $transcripts \
                   --protein_alignments $proteins \
                   --segmentSize 100000 \
                   --overlapSize 10000 \
                   --CPU 32
