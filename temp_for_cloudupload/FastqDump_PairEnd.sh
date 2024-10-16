#!/bin/bash
#SBATCH --time=0-24:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL


module load nixpkgs/16.09

module load sra-toolkit/2.9.6


dir=${1}

for file in $dir*; do
	fastq-dump --split-files $file
done
