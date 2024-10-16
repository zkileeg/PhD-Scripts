#!/bin/bash
#SBATCH --time=0-1:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-44


module load StdEnv/2020  gcc/9.3.0
module load transdecoder/5.7.1

gtf_reference_folder=${1}
fasta_reference_folder=${2}
output_transcript_dir=${3}
#library_name=${4}

log_out=transdecoder_pipe_log_${SLURM_ARRAY_TASK_ID}.txt

echo -e "Input arguments"\n > $log_out
echo -e "gtf_folder: ${gtf_reference_folder}\n" >> $log_out
echo -e "fasta_reference_folder: ${fasta_reference_folder}\n" >> $log_out

echo -e "Getting list of genome files from reference folder\n" >> $log_out
if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
#get list of genome names that have extension .fasta or .fna, both of which are used for genomes
	find $fasta_reference_folder -maxdepth 1 -type f \( -iname \*.fasta -o -iname \*.fna \) | sort > genome_files.txt
	find $gtf_reference_folder -maxdepth 1 -type f -name "*.gtf" | sort > gtf_files.txt
fi


#set genome name dependent on line number from internal array ID number
input_genome=$( awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' genome_files.txt )
input_genome_name="${input_genome##*/}"



#find list of gtf files from input folder
#find $gtf_reference_folder -maxdepth 1 -type f -name "*.gtf" | sort > gtf_files.txt

#get input from list
input_gtf=$( awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' gtf_files.txt )
#input_gtf=$( grep ${input_genome_name} gtf_files.txt )
stringtie_gff="${input_gtf##*/}.gff"


output_transcripts=${output_transcript_dir}${input_genome_name}
output_transcripts_transdecoder="${output_transcripts}.transdecoder.gff3"



echo -e "Running trandecoder with: \n gtf_file=${input_gtf}\ngenome=${input_genome}\n" >> $log_out
echo -e  "Predicting transcripts using transdecoder\n" >> $log_out



#take input as gtf file from stringtie, input fasta genome, and output trancript prediction as fasta
gtf_genome_to_cdna_fasta.pl $input_gtf $input_genome > "${output_transcripts}.fasta"

echo -e "Converting gtf to gff\n" >> $log_out
#take as input gtf file from stringtie and output as gff
gtf_to_alignment_gff3.pl $input_gtf > $stringtie_gff

echo -e "Finding ORFs using the transcripts\n" >> $log_out
#find ORFs using transcripts in fasta format and output as a gff3 file
TransDecoder.LongOrfs -t "${output_transcripts}.fasta"
TransDecoder.Predict -t "${output_transcripts}.fasta"

echo -e "Predicting transcripts using alignment data and ORFs\n" >> $log_out
#get output predicted as gff3 file for the genome
cdna_alignment_orf_to_genome_orf.pl "${input_genome_name}.fasta.transdecoder_dir/longest_orfs.cds.best_candidates.gff3" \
	$stringtie_gff "${output_transcripts}.fasta" > "${input_genome_name}.genome.gff3"
echo -e "Outputting results to ${input_genome_name}.genome.gff3" >> $log_out
