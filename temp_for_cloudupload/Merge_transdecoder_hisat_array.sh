#!/bin/bash
#SBATCH --time=0-00:10:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-44




transdecoder=${1}
hisat=${2}
folder=${3}




