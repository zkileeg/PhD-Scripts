#!/bin/bash
#SBATCH --time=0-06:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL

file_to_compress=${1}

pigz -9 -p32 $file_to_compress
