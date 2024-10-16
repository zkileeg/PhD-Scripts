#!/bin/bash
#SBATCH --time=0-12:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=377G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020  gcc/9.3.0
module load exonerate/2.4.0

#get input from cmd line
protein_reference=${1}
target_assembly=${2}
output_name=${3}

echo $target_assembly

exonerate -m protein2genome -c 1 --querytype protein --targettype dna --percent 70 --minintron 10 --maxintron 60000 --showtargetgff TRUE --showvulgar no --showalignment no \
 $protein_reference $target_assembly > "${output_name}.gff"
