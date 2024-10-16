#!/bin/bash
#SBATCH --time=1-12:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL

module load StdEnv/2023
module load apptainer/1.2.4

genome_in=${1}
cds_in=${2}
WORKDIR=${3}

#run singularity through apptainer. Uses the suggested settings from the EDTA tutorial. 
apptainer run -C -W${SLURM_TMPDIR} -B${WORKDIR}:/home -B /scratch -B /localscratch /scratch/zkileeg/Ecotype_Data/EDTA/EDTA.sif \
	EDTA.pl \
	--genome $genome_in \
	--cds $cds_in \
	--species "others" \
	--step "all" \
	--overwrite 1 \
	--sensitive 1 \
	--anno 1 \
	--evaluate 1 \
	--threads 32

#there are a bajillion files. As file numbers are capped, remove the intermediate files we don't need. 
rm -r "$WORKDIR"zkileeg/*raw
rm -r "$WORKDIR"zkileeg/*final
rm -r "$WORKDIR"zkileeg/*combine
