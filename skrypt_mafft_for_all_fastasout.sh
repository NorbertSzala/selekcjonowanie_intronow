#!/bin/bash

#opis:
#skrypt robiacy maffta na zadanych parametrach z outputem w formacie fasta w danym folderze


input_path="/home/norbert/mrdn/euglena/kod_i_pliki/surowe_pliki_plus_minus_500/raw_reads_9/merging_fastas"
output_path="/home/norbert/mrdn/euglena/kod_i_pliki/surowe_pliki_plus_minus_500/raw_reads_9/after_mafft_fastas"

# Sprawdzenie czy folder wejściowy istnieje
if [ ! -d "$input_path" ]; then
    echo "Folder wejściowy $input_path doesn't exist"
    exit 1
fi

# Sprawdzenie czy folder wyjściowy istnieje
if [ ! -d "$output_path" ]; then
    echo "Folder $output_path doesn't exist. Creating folder..."
    mkdir -p "$output_path" || { echo "Error. Creating folder didn't succeed."; exit 1; }
fi



# Przetwarzanie plików fasta
for fasta_file in "$input_path"/*fasta; do
    if [ -f "$fasta_file" ] && [ -s "$fasta_file" ]; then
        filename=$(basename "$fasta_file")
        filename_without_extension="${filename%.*}"  # Usuń rozszerzenie .fasta

        # Usuwanie znaków "-" tylko z linii niezaczynających się od ">"
        while IFS= read -r line; do
            if [[ $line != ">"* ]]; then
                cleaned_line=$(echo "$line" | tr -d '-')
                echo "$cleaned_line" >> "$output_path/$filename_without_extension.tmp"
            else
                echo "$line" >> "$output_path/$filename_without_extension.tmp"
            fi
        done < "$fasta_file"

        # Wykonanie mafft na przetworzonym pliku
        mafft --op 1.53 --ep 0.5 --jtt 1 --thread 50 "$output_path/$filename_without_extension.tmp" > "$output_path/$filename_without_extension.aln"

         #Usunięcie pliku tymczasowego
        rm "$output_path/$filename_without_extension.tmp"
    else
        echo "File $fasta_file doesn't exist or is empty"
    fi
done