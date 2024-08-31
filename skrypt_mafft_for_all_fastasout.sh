#!/bin/bash

#opis:
#skrypt robiacy maffta na zadanych parametrach z outputem w formacie fasta w danym folderze


input_path='/home/norbert/mrdn/euglena/kod_i_pliki/surowe_pliki_plus_minus_500/raw_reads_9/rozdzielone_pliki/'
output_path='/home/norbert/mrdn/euglena/kod_i_pliki/surowe_pliki_plus_minus_500/raw_reads_9/pliki_po_maffcie/'

#input_path="/home/szala/euglena/kod_i_pliki/divided_fastas"
#output_path="/home/szala/euglena/kod_i_pliki/fastas_after_mafft"

if [ ! -d "$input_path" ]; then
    echo "Folder wejściowy $input_path nie istnieje"
    exit 1
fi

# Sprawdzenie, czy folder wyjściowy istnieje
if [ ! -d "$output_path" ]; then
    echo "Folder $output_path nie istnieje. Tworzenie folderu..."
    mkdir -p "$output_path" || { echo "Błąd. Tworzenie folderu nie powiodło się."; exit 1; }
fi

# Przetwarzanie plików FASTA
for fasta_file in "$input_path"/*.fasta; do
    if [ -f "$fasta_file" ] && [ -s "$fasta_file" ]; then
        filename=$(basename "$fasta_file")

        tmp_file="$output_path/filename.tmp"

        # Usuwanie znaków "-" tylko z linii niezaczynających się od ">"
        while IFS= read -r line; do
            if [[ $line != ">"* ]]; then
                cleaned_line=$(echo "$line" | tr -d '-')
                echo "$cleaned_line" >> "$tmp_file"
            else
                echo "$line" >> "$tmp_file"
            fi
        done < "$fasta_file"

        # Wykonanie MAFFT na przetworzonym pliku
        mafft --op 1.53 --ep 0.5 --jtt 1 --thread 50 "$tmp_file" > "$output_path/filename.fasta"

        # Usunięcie pliku tymczasowego
        rm "$tmp_file"
    else
        echo "Plik $fasta_file nie istnieje lub jest pusty"
    fi
done
