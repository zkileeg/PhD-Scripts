#!/bin/bash
#SBATCH --time=0-04:0:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL

#load necessary modules
module load StdEnv/2020
module load minimap2/2.24
module load racon/1.4.13

#get positional input from command line. Input and output are file names, no extensions.
fasta_assembly=${1}
fastq_in=${2}
output=${3}
sam_out=${4}

#iterate three times, mapping the long-reads to our genome, polishing with racon, then doing this again on the newly polished genome each time
flag=1
while [ $flag -le 3 ]; do

	sam_assembly="${sam_out}${flag}.sam"
	racon_output="${output}${flag}.fasta"

	minimap2 -t 32 -ax map-ont $fasta_assembly $fastq_in > $sam_assembly



	racon -t 32 $fastq_in $sam_assembly $fasta_assembly > $racon_output

	fasta_assembly=$racon_output



	flag=$(( $flag + 1 ))


done
