#!/bin/bash


input_path="/home/norbert/mrdn/euglena/kod_i_pliki/surowe_pliki_plus_minus_500/raw_reads_9/merging_fastas"
output_path="/home/norbert/mrdn/euglena/kod_i_pliki/surowe_pliki_plus_minus_500/raw_reads_9/after_mafft_fastas"

for fasta_file in "$input_path"/*fasta; do
	if [ -f "$fasta_file" ] && [ -s "$fasta_file" ]; then
		filename=$(basename "$fasta_file")
		filename_without_extension="${filename%.*}"  # Usuń rozszerzenie .fasta	
		mafft --op 1.53 --ep 0.5 --jtt 1 --clustalout  "$fasta_file" > "$output_path/$filename_without_extension.aln"
	else
		echo "File $fasta_file doesn't exist or is empty"
	fi
done
