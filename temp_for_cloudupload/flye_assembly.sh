#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL


#run flye assembly with three iterations of polishing

#load necessary modules
load module gcc python/3.1.0

#load python virtual environment so we can run flye
load virtual environment
#load environment from where it's saved 
source flye/bin/activate

input_reads=${1}
genome_size=${2}
outdir=${3}

#run flye under raw mode, meaning flye will do some polishing steps 
python ~/Flye-2.9.2/bin/flye --pacbio-raw $input_reads  --threads 32 --iterations 3 -o $outdir --genome-size $genome_size
