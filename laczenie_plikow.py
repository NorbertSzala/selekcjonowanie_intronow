#!/usr/bin/env python3

import os

#serwer
path_input = '/home/milewski/eugleny/Mafft2_addlong'
path_output = '/home/szala/euglena/kod_i_pliki/divided_fastas'

#lokalnie
# path_input = '/home/norbert/mrdn/euglena/kod_i_pliki/surowe_pliki_plus_minus_500/raw_reads_9/zlaczone_pliki'
# path_output = '/home/norbert/mrdn/euglena/kod_i_pliki/surowe_pliki_plus_minus_500/raw_reads_9/rozdzielone_pliki'

# merging fasta files into one file
#If you would like to overwrite a files, first delete them manually
#uporzadkujmy tak by jednak ta sekwencja z wieksza iloscia przerw byla pierwsza #python nie obsluguje kolejnosci w
# slowniku - bedziemy musieli to uwzglednic w kodzie na funkcje cutting_scrap # we wszystkich plikach jakie recznie
# sprawdziłem tak własnie jest i nie jest wymagana modyfikacja

def merging_fastas(path_input, path_output):
    for filename in os.listdir(path_input):
        file = os.path.join(path_input, filename)
        choosen_organism = None
        sequences_dict = {}
        if file.endswith(".fasta"):
            with open(file, "r") as file:
                for line in file:
                    if line.startswith(">"):
                        organism_id = line.strip().replace(" ", "_") #that first fine in fasta file
                        choosen_organism = organism_id
                        if choosen_organism not in sequences_dict:
                            sequences_dict[choosen_organism] = []
                    else:
                        sequences_dict[choosen_organism].append(line.strip())

            for name_from_dictionary in sequences_dict.keys():
                print("name from dictionary:", name_from_dictionary)
                print("filename: ", filename)
                output_file = f"{name_from_dictionary[1:4]}_{str(filename)}"
                path_to_output_file = os.path.join(path_output, output_file)
                print("output: ", output_file, "\n")
                if not os.path.exists(path_to_output_file):
                    print(f"tworze ten plik {path_to_output_file} poraz pierwszy")
                    with open(path_to_output_file, "w") as f:
                        f.write(name_from_dictionary + "\n")
                        for lines in sequences_dict[name_from_dictionary]:
                            f.write(lines + "\n")
                else:
                    count = 0
                    with open(path_to_output_file, "r") as f:
                        text = f.read()
                        count = text.count(">")
                        if count < 2:
                            with open(path_to_output_file, "a") as f:
                                f.write(name_from_dictionary + "\n")
                                for lines in sequences_dict[name_from_dictionary]:
                                    f.write(lines + "\n")
                        else:
                            print(f"In that file {path_to_output_file} there are already 2 fasta sequences. "
                                  f"Make sure is it correct")


merging_fastas(path_input, path_output)