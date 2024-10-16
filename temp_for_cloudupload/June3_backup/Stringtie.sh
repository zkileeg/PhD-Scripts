#!/bin/bash
#SBATCH --time=0-12:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020
module load stringtie/2.2.1

bam_in=${1}
gtf_output=${2}

stringtie -o $gtf_output -p 32 $bam_in
