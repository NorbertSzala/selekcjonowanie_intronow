#!/usr/bin/env python3


from Bio import AlignIO, Align
import pandas as pd
import os
import numpy as np
import time

#description of goal
    #main goal that function is select introns and exon from alignment and descripte their class, start - ends nucleotide.
#description step by step
    # 1. iterate through files in input before mafft and linking with input after mafft, by the same name
    # 2. checking if both files have .fasta extension
    # 3. making alignment in file after mafft
    # 4. extracting_strands_from_alignment function is used to determine which strand is genomic strand or transcript
    # 5. Removing gaps from 5' and 3'. Then pre-exon (candidat to be exon) will be first
    # 6. Creating list with nucleotides's positions of intron or exon
    # 7. Creating dictionary like: start_exon:end_exon. If distance between two following indices is bigger than two (acceptable gap length) then there is intron between them.
    # 8. Making alignment and determinig class of exon according to their homology.
    # 9. Forming a table with exons
    # 10. Spliting a data frame to 3 others, with different exon's class
    # 11. Adding introns to tables
    # 12. Saving table in tsv gff format.

#ASSUMPTIONS
#- both files, path_to_file_before_MAFFT and path_to_file_after_MAFFT have to have the same name
#- they have to be .fasta
#- acceptable_gap_length is int.

min_length_aligned_sequence = 30 #Minimal lenght of sequence which could be an exon
# extreme_homology = 0.97 #percentage of homology of sequence, treshold #I assume two faulty aligned nucleotides per 100 (98%) and one more nt because sometimes latest nt can move from end of one sequence to beginning next sequence

path_to_file_before_MAFFT = "/home/szala/euglena/kod_i_pliki/divided_fastas"
path_to_file_after_MAFFT = "/home/szala/euglena/kod_i_pliki/fastas_after_mafft"

def cutting_scrap(path_to_file_before_MAFFT, path_to_file_after_MAFFT, acceptable_gap_length, extreme_homology):
    gff_data_frames = [] #pre-list for tables
    gaps_signs = "-" * acceptable_gap_length #maximum length of gaps in sequence in one exon's sequence
    files_in_progress = 0
    broken_files = [] #list of files which caused error
    nucleotides = ["a", "t", "g", "c", "A", "T", "G", "C"]
    column_names = ['seqid', 'source', 'exon_type', 'start', 'end', 'length', 'homology', 'strand', 'phase', 'attributes']
    print(f"\n \n \n \n Start running function: cutting_scrap")

    total_files = len([f for f in os.listdir(path_to_file_after_MAFFT) if f.endswith(".fasta")])

    #linking two files from input and making alignment
    for filename in os.listdir(path_to_file_after_MAFFT):
        file = os.path.join(path_to_file_after_MAFFT, filename)
        if not os.path.isfile(file):
            continue
        if file.endswith(".fasta"):
            try:
                files_in_progress += 1
                count_files = percentage_of_advancement(path_to_file_before_MAFFT)
                print(f"\n ################################################################################### \n running file: {filename} which is {files_in_progress} of {count_files}, what means {round(files_in_progress/count_files*100, 2)}% of advancement \n ################################################################################### \n")
                alignment = AlignIO.read(file, "fasta")
                alignment_id = alignment[0].id
                alignment_DIDNT_TOUCHED = AlignIO.read(file, "fasta")
    
    
                file_before_MAFFT = os.path.join(path_to_file_before_MAFFT, filename)
                if not os.path.isfile(file_before_MAFFT):
                    continue
        
                # extracting_strands_from_alignment function is used to determine which strand is genomic strand or transcript
                seq1_DT, seq2_DT = extracting_strands_from_alignment(alignment_DIDNT_TOUCHED)
                seq1, seq2 = extracting_strands_from_alignment(alignment)
    
                # Removing gaps from 5' and 3'. Then pre-exon (candidat to be exon) will be first
                seq1, seq2, count_deleted_gaps_left = cleaning_gaps_from_both_edges(seq1_DT, seq2_DT, alignment, filename) #deleting gaps from 5' and 3'. It is important. The next function start counting from possible exon.
        
                # Creating list with nucleotides's positions of intron or exon
                temp_intron, temp_exon, invalid_nucleotides = indices_of_introns_and_exons_fun(seq1_DT, seq2_DT, seq1, seq2, count_deleted_gaps_left, alignment_id)
                
                # Creating dictionary like: start_exon:end_exon. If distance between two following indices is bigger than two (acceptable gap length) then there is intron between them.
                all_exon_range_dict = start_and_end_parameters_of_exons_dict_fun(temp_exon, gaps_signs)
        
        
                # Making alignment and determinig class of exon according to their homology.
                # Forming a table with exons
                exon_rows_to_concete_table, one_alignment_df = making_exons_gff_table_fun(all_exon_range_dict, min_length_aligned_sequence, extreme_homology, seq1_DT, seq2_DT, alignment_id, column_names, filename[:-6])

                if one_alignment_df.empty:
                    print(f'File: {filename} is out of exons and introns. Processing failed. Broken_file list appended.')
                    broken_files.append(filename)
                else:
                    gff_data_frames.append(one_alignment_df)
            
            except Exception as e:
                print(f"Error processing file {filename}: {e}")
                broken_files.append(filename)
        
    
    # spliting a data frame to 3 others, with different exon's class
    if not len(gff_data_frames) == 0:
        all_exons_df = pd.concat(gff_data_frames, ignore_index = True)
        TLH_exons_df = all_exons_df[all_exons_df['exon_type'] == 'TLH'].copy()
        fine_exons_df = all_exons_df[all_exons_df['exon_type'] == 'fine'].copy()

        all_exons_df, TLH_exons_df, fine_exons_df = adding_introns_to_gff_data_frame(all_exons_df, TLH_exons_df, fine_exons_df, column_names)
  
        print(f'{((total_files - len(broken_files)) / total_files) * 100 }% of files end up succesfully. {len(broken_files)} had error. List of files with errors in broken_files file. Invalid nucleotides (different than A T G C) are written with positions in invalid_nucleotides file')
        fine_exons_count = fine_exons_df['exon_type'].str.contains('fine', case = False).sum()

    else:
        print("No valid data frames to concatenate.")
        all_exons_df = pd.DataFrame(columns=column_names)
        TLH_exons_df = pd.DataFrame(columns=column_names)
        fine_exons_df = pd.DataFrame(columns=column_names)
        fine_exons_count = 0        
        
    # Adding introns to tables

    # Saving table in tsv gff format
    # save_to_gff_file(all_exons_df, 'broken_all_exons_gff.tsv')
    # save_to_gff_file(TLH_exons_df, 'broken_TLH_exons_gff.tsv')
    # save_to_gff_file(fine_exons_df, 'broken_fine_exons_gff.tsv')
    print(len(broken_files))
    with open('broken_files', 'w') as file:
        file.write('lista plikow ktora nie przeszla programu \n ')
        file.write(str(broken_files))
    length_fine_exons_df = len(fine_exons_df)

    with open('invalid_nucleotides', 'w') as file:
        for alignment_id, nt1, nt2, pos in invalid_nucleotides:
            file.write(f'{alignment_id}\t{nt1}\t{nt2}\t{pos}\n')
        
    return fine_exons_count, length_fine_exons_df

                   

#############################################################################################################################
#######################################                MINOR FUNCTIONS                #######################################       
#############################################################################################################################

def percentage_of_advancement(directory):
    count = 0
    for file in os.listdir(directory):
        if file.endswith(".fasta"):
            count += 1
    return(count)
    

def cleaning_gaps_from_both_edges(sequence_transcipt, sequence_genome, alignment, filename):
    sequence_transcipt_left_shorted = sequence_transcipt.lstrip("-")
    count_deleted_gaps_left = len(sequence_transcipt) - len(sequence_transcipt_left_shorted) #inform how many gaps were deleted from 5'
    
    sequence_transcipt_right_and_left_shorted = sequence_transcipt_left_shorted.rstrip("-")
    count_deleted_gaps_right = len(sequence_transcipt_left_shorted) - len(sequence_transcipt_right_and_left_shorted) #inform how many gaps were deleted from 3'
    
    if count_deleted_gaps_right:
        sequence_genome_left_and_right_shorted = sequence_genome[count_deleted_gaps_left:-count_deleted_gaps_right]
    else:
        sequence_genome_left_and_right_shorted = sequence_genome[count_deleted_gaps_left:]

    if len(sequence_genome_left_and_right_shorted) != len(sequence_transcipt_right_and_left_shorted):
        raise ValueError(f" VALUE ERROR: Given sequences: {alignment[0].id} - len: {len(alignment[0])} and  {alignment[1].id} - len {len(alignment[1])} from {filename} must have the same length.")
    return sequence_transcipt_right_and_left_shorted, sequence_genome_left_and_right_shorted, count_deleted_gaps_left


def indices_of_introns_and_exons_fun(seq1_DT, seq2_DT, seq1, seq2, count_deleted_gaps_left, alignment_id):
    #linking pre-exon's nucleotides (nt-nt pairs) and pre-intron's nucleotides (gap-nt pairs) to theirs indices 
    #NOTE: output is nucleotide or gap from seq1 (transcript sequence). Even if it is gap like this: 
    #seq1: aaatttggg, seq2: aaa---ggg, output will be 'ttt' instead '---'
    #it shows that where were gap, seq2 genome or in seq1 transcript
    temp_exon = [] #list containing indices of intron's positions in sequence.
    temp_intron = [] #list containing indices of intron's positions in sequence.
    invalid_nucleotides = []
    index = 0
    paired_nucleotides = zip(seq1, seq2)
    valid_nucleotides = set('ATGCatgc')
    
    if len(seq1_DT) != len(seq2_DT):
        raise ValueError("The two sequences must be of the same length.")
        print(seq1)
        print(seq2)
        
    for nt1, nt2 in paired_nucleotides:
        if "-" in (nt1, nt2):
            temp_intron.append(index + count_deleted_gaps_left)
        elif nt1 not in valid_nucleotides or nt2 not in valid_nucleotides:
            if alignment_id not in invalid_nucleotides:
                invalid_nucleotides.append((alignment_id, nt1, nt2, index + count_deleted_gaps_left))  # Dodanie nukleotydów i pozycji
        else:
            temp_exon.append(index + count_deleted_gaps_left)
        index += 1 #indeks do wskazywania pozycji w sekwencji 
    return temp_intron, temp_exon, invalid_nucleotides


def start_and_end_parameters_of_exons_dict_fun(temp_exon, gaps_signs):
    #linking indices into single strand of possible exons
    start_exon = temp_exon[1]
    all_exon_range_dict = {} #that dictionary contains: key(index of start exon) and value(index of end exon)
    
    for i in range(len(temp_exon)): 
        #print(start_exon)
        if (temp_exon[i] - temp_exon[i-1]) > 1+len(gaps_signs): #if difference between two indices of exons's positions is bigger than given number (2), that smaller number is index of 3' nucleotide in exon 
            end_exon = temp_exon[i-1] #temp_exon[i-1] because it is index in list temp_exon. +1 because python starts counting from 0
            all_exon_range_dict[start_exon] = end_exon +1  #creating dictionary with all exons, even with theese too short and with too low homology
            start_exon = temp_exon[i] +1
            #brakuje ostatniego eksonu. Jest to spowodowane tym, ze nie mozna okreslic parametru end_exon, gdyz nie w liscie temp_exon wiekszego numeru niz ten ostatni
    all_exon_range_dict[start_exon] = temp_exon[-1] #last pair

    if all_exon_range_dict == {}:
        raise ValueError(f"ValueError - all_exon_range_dict seems to be empty. Its length = {len(all_exon_range_dict)}")
    #print(f" \n all_exon_range_dict: {all_exon_range_dict}")
    return (all_exon_range_dict)


     
def start_and_end_parameters_of_introns_dict_fun(all_exon_range_dict):
    #linking indices into single strand of possible introns
    intron_range_dict = {}
    keys_from_all_exon_range_dict = sorted(all_exon_range_dict.keys()) #list that contains exons's start positions
    #print(f"keys_from_all_exon_range_dict: {keys_from_all_exon_range_dict}\n")

    for i in range(len(keys_from_all_exon_range_dict) - 1):
        start_intron = all_exon_range_dict[keys_from_all_exon_range_dict[i]] +1 
        end_intron = keys_from_all_exon_range_dict[i + 1] - 1
        intron_range_dict[start_intron] = end_intron
    if intron_range_dict == {}:
        raise ValueError(f"ValueError - intron_range_dict seems to be empty. Its length = {len(intron_range_dict)}")
    #print(f" \n intron_range_dict: {intron_range_dict}")
    return (intron_range_dict)


#zmienione all_fine_exons na fine_exons
def making_exons_gff_table_fun(all_exon_range_dict, min_length_aligned_sequence, extreme_homology, seq1_DT, seq2_DT, alignment_id, column_names, filename):
    exon_rows_to_concete_table = []

    aligner = Align.PairwiseAligner()
    aligner.mismatch_score = 0 #customized setting towards get pure percentage of homology, not alignment with gap penalty etc. The object of interest is simply that how many nt in query has equivalend in target seq. 
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0
    
    for key in all_exon_range_dict:    
        query = seq1_DT[key -1 : all_exon_range_dict[key]]
        target = seq2_DT[key-1 : all_exon_range_dict[key]]
        target_length = len(target)
    
        if target_length > min_length_aligned_sequence:
            score = aligner.score(query, target)
            identity_level = (round(((score * 100) / len(target)), 2))
            #print(f"query: {query}, \n target: {target}, percent of homology = {round(score / len(target), 2) * 100} %, score = {score}, len = {len(target)}\n")
            if identity_level >= extreme_homology * 100:
                exon_class = 'fine'
                # print(f"\n query: {query} \n target: {target} \n with start: {key} and end: {all_exon_range_dict[key]} goes to fine_exons with homology: {identity_level} \n")
                # if len(target) != len(query):
                    # print(f'target: {len(target)} query: {len(query)}')
            else:
                exon_class = 'TLH'  
                # print(f"\n query: {query} \n target: {target} \n with start: {key} and end: {all_exon_range_dict[key]} goes to TLH_exons with homology: {identity_level}")
                # if len(target) != len(query):
                    # print(f'target: {len(target)} query: {len(query)}')
        else:
            identity_level = 0
            exon_class = "TS" 
            # print(f"\n query: {query} \n target: {target} \n with start: {key} and end: {all_exon_range_dict[key]} goes to TS_exons with homology: {identity_level}")
            # if len(target) != len(query):
                    # print(f'target: {len(target)} query: {len(query)}')
        source = extending_source_data_frame(alignment_id)
        strand = extending_strand_data_frame(alignment_id)
        exon_rows_to_concete_table.append((alignment_id, source, exon_class, (key), (all_exon_range_dict[key]), ((all_exon_range_dict[key]) - (key)+1), identity_level, strand, ".", filename))
    df = pd.DataFrame(exon_rows_to_concete_table, columns = column_names)
    
    return exon_rows_to_concete_table, df

    

def extending_source_data_frame(alignment_id):
    #print(f"Alignment: {alignment[0].id}, funkcja: extending_source_data_frame, przeszlo")
    keywords = ["TRINITY", "BACKBONE", "SCAFFOLD"]
    found = False
    for key in keywords:
        alignment_id = str(alignment_id).lower()
        if alignment_id.find(key.lower()) != -1:
            source = key
            found = True
            break
    if not found:
        source = None #none means unknown
    return source


def extending_strand_data_frame(alignment_id): 
    #print(f"Alignment: {alignment[0].id}, funkcja: extending_strand_data_frame, przeszlo")
    if str(alignment_id).find("SL+") != -1:
        strand = "+"
    elif str(alignment_id).find("SL-") != -1:
        strand = "-"
    else:
        strand = None
    return strand


def adding_introns_to_gff_data_frame(all_exons_df, TLH_exons_df, fine_exons_df, column_names):
    data_frames_list = [all_exons_df, TLH_exons_df, fine_exons_df]

    for i, df in enumerate(data_frames_list):
        df.sort_values(by = ['seqid', 'start'], inplace = True)
        
        df['seqid_introns'] = df['seqid'].shift(-1) #creating new column with seqid of next sequention. In next steps, rows without matching seqid and seqid_introns, would not be consideresd
        df['intron_start'] = df['start'].shift(-1) #creating index of start position in introns
        
        intron_df_temp = df[df['seqid'] == df['seqid_introns']].copy()
        intron_df_temp['start'] = df['end'] + 1
        intron_df_temp['end'] = df['intron_start'] - 1
        intron_df_temp['exon_type'] = 'intron'
        intron_df_temp['homology'] = 0
        intron_df_temp['length'] = intron_df_temp['end'] - intron_df_temp['start'] + 1
        intron_df_temp = intron_df_temp[column_names]

        df = df[column_names]        
        df = pd.concat([df, intron_df_temp]).sort_values(by = ['seqid', 'start']).reset_index(drop = True)

        data_frames_list[i] = df

    return data_frames_list[0], data_frames_list[1], data_frames_list[2]


def save_to_gff_file(gff_final_data_frame, filename):
    #print(f"Alignment: {alignment[0].id}, funkcja: save_to_gff_file, przeszlo")
    if os.path.isfile(filename):
        user_input = input(f"GFF file  already exists. Do you want overwrite? y/n: ")
        if user_input == "n":
            base_name, ext = os.path.splitext(filename)
            i = 1
            while os.path.isfile(f"{base_name}_{i}{ext}"):
                i += 1
            filename = f"{base_name}_{i}{ext}"
            gff_final_data_frame.to_csv(filename, sep = "\t")
            
        elif user_input == "y":
            gff_final_data_frame.to_csv(filename, sep = "\t")
        else:
            print("Clarify your answer. Nothing has done.")
            
    else:
        gff_final_data_frame.to_csv(filename, sep = "\t")


def extracting_strands_from_alignment(alignment): #that function describe what strand is transcript strand and genomic strand - in our files, transcript strand has longer ID.
    #print(f"Alignment: {alignment[0].id}, funkcja: extracting_strands_from_alignment, przeszlo")
    seq1 = alignment[0].seq
    seq2 = alignment[1].seq
    # if len(alignment[0].id) >= len(alignment[1].id):
    #     seq1 = alignment[0].seq
    #     seq2 = alignment[1].seq
    #     #seq1_is_transcript = True #obecnie nie uzywane
    # else:
    #     seq1 = alignment[1].seq
    #     seq2 = alignment[0].seq    
        #seq1_is_transcript = False #obecnie nie uzywane
    #print(f"\n Just run function: extracting_strands_from_alignment")
    return seq1, seq2
    #return seq1_is_transcript #obecnie nie uzywane

def function_progress(alignment, alignment_DIDNT_TOUCHED, files_in_progress, count_files):
    #do tego tematu wrocimy jak zrobimy tabele z odrzuconymi wartosciami
    first_50_nt_alignment = str((alignment[1].seq)[:20])
    #print(f" {round((100-(len(alignment[0])\/len(alignment_DIDNT_TOUCHED[0]))*100), 2)}% % base pair done of {alignment[0].id}. {files_in_progress} of {count_files}, what means {round(files_in_progress/count_files*100, 2)}% of advancement. The time is: {time.strftime('%H:%M:%S', time.localtime())}")

if __name__ == '__main__':
    cutting_scrap(path_to_file_before_MAFFT, path_to_file_after_MAFFT, 2, 0.97)