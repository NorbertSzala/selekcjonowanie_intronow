#!/usr/bin/env python3

from Bio import AlignIO, Align
import pandas as pd
import os
import numpy as np
from datetime import datetime
from IPython.display import display, HTML

# IMPORTANT: fine_exons_df is not perfect, better one is final_fine_exon_sub_dfs. In output, there will be only
# final_fine_exon_sub_dfs and all_exons_df

'''
# description of goal main goal that function is select introns and exon from alignment and describe their class,
# start - ends nucleotide. description step by step:
# 1. iterate through files in input before mafft and linking with input after mafft, by the same name
# 2. checking if both files have .fasta extension
# 3. Extracting_strands_from_alignment
# 4. Removing gaps from 5' and 3'. Then pre-exon (candidate to be exon) will be first
# 5. Creating list with nucleotides's positions of intron or exon
# 6. Creating dictionary like: start_exon:end_exon. If distance between two following indices is bigger than two
# (acceptable gap length) then there is intron between them.
# 7. Creating substring in format: exon intron exon XXXX - if XXXX is bad exon, substring equals exon intron exon.
# Another substring is creating by next proper exon. Substring are used to make final_fine_exon_intron_sub_dfs
# 8. Making alignment and determining class of exon according to their homology. It is used to make all_exons_df
# 9. Forming a table with exons
# 10. Splitting a data frame to 3 others, with different exon's class
# 11. Adding introns to tables
# 12. Saving table in tsv gff format.
'''

# ASSUMPTIONS
# - both files, path_to_file_before_mafft and path_to_file_after_mafft have to have the same name
# - they have to be .fasta
# - acceptable_gap_length is int.

# min_length_aligned_sequence = 30  # Minimal length of sequence which could be an exon
# extreme_homology = 0.97 #percentage of homology of sequence, threshold

# na serwer
path_to_file_before_mafft = "/home/szala/euglena/kod_i_pliki/divided_fastas"
path_to_file_after_mafft = "/home/szala/euglena/kod_i_pliki/fastas_after_mafft"



# lokalnie
# path_to_file_before_mafft = ('/home/norbert/mrdn/euglena/kod_i_pliki/surowe_pliki_plus_minus_500/raw_reads_9/'
#                              'merging_fastas')
# path_to_file_after_mafft = ('/home/norbert/mrdn/euglena/kod_i_pliki/surowe_pliki_plus_minus_500/raw_reads_9/'
#                             'fastas_after_mafft_na_probe')


#lokalnie caly zbior
# path_to_file_before_mafft = '/home/norbert/mrdn/euglena/kod_i_pliki/Nowa_wersja_funkcji/divided_fastas'
# path_to_file_after_mafft = '/home/norbert/mrdn/euglena/kod_i_pliki/Nowa_wersja_funkcji/fastas_after_mafft'


def cutting_scrap(path_to_file_before_mafft, path_to_file_after_mafft, acceptable_gap_length, extreme_homology,
                  min_length_aligned_sequence):
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', None)

    now = datetime.now()
    current_time = now.strftime("%Y-%m-%d_%H:%M")
    gff_data_frames = []  # pre-list for tables
    final_fine_exon_sub_df = []
    final_fine_intron_sub_df = []
    invalid_nucleotides = []
    average_introns_length = {}
    unique_seqid_in_final_df_number_all = {}
    substrings_introns_number_all = {}
    files_in_final_fine_exon_intron_sub_dfs_all = {}
    introns_count_all = {}
    exons_count_all = {}

    gaps_signs = "-" * acceptable_gap_length  # maximum length of gaps in sequence in one exon's sequence
    files_in_progress = 0
    alignments_with_one_exon = 0
    broken_files = []  # list of files which caused error
    column_names = ['seqid', 'source', 'exon_type', 'start', 'end', 'length', 'homology', 'strand', 'phase',
                    'attributes']
    column_names_substring = ['seqid', 'source', 'exon_type', 'start', 'end', 'length', 'homology', 'GC', 'strand', 'phase',
                              'attributes', 'substring']
    species = ['gra', 'lon', 'hie']
    print(f"\n \n \n \n Start running function: cutting_scrap")

    total_files = len([f for f in os.listdir(path_to_file_after_mafft) if f.endswith(".fasta")])

    # linking two files from input and making alignment
    for filename in os.listdir(path_to_file_after_mafft):
        file = os.path.join(path_to_file_after_mafft, filename)
        if not os.path.isfile(file):
            continue
        if file.endswith('.fasta'):
            try:
                files_in_progress += 1
                count_files = percentage_of_advancement(path_to_file_before_mafft)
                print(f"\n ################################################################################### \n "
                      f"running file: {filename} which is {files_in_progress} of {count_files}, what means "
                      f"{round(files_in_progress / count_files * 100, 2)}% of advancement \n "
                      f"################################################################################### \n")
                alignment = AlignIO.read(file, "fasta")
                alignment_id = alignment[0].id
                alignment_didnt_touched = AlignIO.read(file, "fasta")

                file_before_mafft = os.path.join(path_to_file_before_mafft, filename)
                if not os.path.isfile(file_before_mafft):
                    continue

                # extracting_strands_from_alignment function is used to determine which strand is genomic strand
                # or transcript
                seq1_dt, seq2_dt = extracting_strands_from_alignment(alignment_didnt_touched)
                seq1, seq2 = extracting_strands_from_alignment(alignment)

                # Removing gaps from 5' and 3'. Then pre-exon (candidat to be exon) will be first

                seq1, seq2, count_deleted_gaps_left = cleaning_gaps_from_both_edges(seq1_dt, seq2_dt, alignment,
                                                                                    filename)  # deleting gaps from 5'
                # and 3'. It is important. The next function start counting from possible exon.

                # Creating list with nucleotides's positions of intron or exon
                temp_intron, temp_exon, invalid_nucleotides = indices_of_introns_and_exons_fun(seq1_dt, seq2_dt, seq1,
                                                                                               seq2,
                                                                                               count_deleted_gaps_left,
                                                                                               alignment_id,
                                                                                               invalid_nucleotides)

                # Creating dictionary like: start_exon:end_exon. If distance between two following indices is bigger
                # than two (acceptable gap length) then there is intron between them.
                all_exon_range_dict = start_and_end_parameters_of_exons_dict_fun(temp_exon, gaps_signs)

                # method with higher accuracy. Introns exist when they are surrounded by two exons like: ex int ex.
                # Create substrings with only fine exons. Sometimes in one alignment there it more than one substring,
                # when between two fine exons is one bad exon
                fine_exon_sub_dict, fine_intron_sub_dict, alignments_with_one_exon = create_substring_dict(
                    all_exon_range_dict, seq1_dt, seq2_dt, min_length_aligned_sequence, extreme_homology,
                    alignments_with_one_exon)

                # creating data frame with fine exons
                fine_exon_sub_df = create_final_exon_gff_table_from_substring(fine_exon_sub_dict, seq1_dt, seq2_dt,
                                                                              alignment_id, column_names_substring,
                                                                              filename[:-6])
                fine_intron_sub_df = create_final_intron_table_from_substring(fine_intron_sub_dict, alignment_id,
                                                                              column_names_substring, filename[:-6], seq2_dt)

                # Making alignment and determining class of exon according to their homology.
                # Forming a table with exons
                exon_rows_to_concat_table, one_alignment_df = making_exons_gff_table_fun(all_exon_range_dict,
                                                                                         min_length_aligned_sequence,
                                                                                         extreme_homology, seq1_dt,
                                                                                         seq2_dt, alignment_id,
                                                                                         column_names, filename[:-6])

                if fine_exon_sub_df.empty:
                    print(f'File: {filename} is out of exons and introns. Processing finished. No_exons_files file '
                          f'extended.')
                else:
                    final_fine_exon_sub_df.append(fine_exon_sub_df)
                    final_fine_intron_sub_df.append(fine_intron_sub_df)

                gff_data_frames.append(
                    one_alignment_df)  # gff data frames with inaccurate introns doesn't report errors

            except Exception as e:
                print(f"Error processing file {filename}: {e}")
                broken_files.append(filename)

    # spliting a data frame to 3 others, with different exon's class

    if not len(final_fine_exon_sub_df) == 0:
        final_fine_exon_sub_dfs = pd.concat(final_fine_exon_sub_df, ignore_index=True)
        final_fine_intron_sub_dfs = pd.concat(final_fine_intron_sub_df, ignore_index=True)

        # final_fine_exon_intron_sub_dfs = adding_introns_to_final_gff_data_frame(final_fine_exon_sub_dfs,
        #                                                                         column_names_substring)
        final_fine_exon_intron_sub_dfs = pd.concat([final_fine_exon_sub_dfs, final_fine_intron_sub_dfs],
                                                   ignore_index=True)

        final_fine_exon_intron_sub_dfs = final_fine_exon_intron_sub_dfs[~((final_fine_exon_intron_sub_dfs['length'] < 40)
                                                                        & (final_fine_exon_intron_sub_dfs['exon_type']
                                                                           == 'intron'))]

        final_fine_exon_intron_sub_dfs = final_fine_exon_intron_sub_dfs.sort_values(by=['attributes', 'start'])

        files_in_final_fine_exon_intron_sub_dfs_all, files_in_final_fine_exon_intron_sub_dfs, unique_seqid_in_final_df_number = (
            getting_files_from_df(final_fine_exon_intron_sub_dfs, 'attributes',
                                  files_in_final_fine_exon_intron_sub_dfs_all, species))



        # Extracing information from final df
        filtered_introns_df = final_fine_exon_intron_sub_dfs[
            final_fine_exon_intron_sub_dfs['exon_type'] == 'intron'].copy()

        introns_count_all, average_introns_length, substrings_introns_number_all, unique_seqid_in_final_df_number_all = extracting_informations_to_return(
            filtered_introns_df, species, introns_count_all, average_introns_length, substrings_introns_number_all,
            unique_seqid_in_final_df_number_all)

        for specie in species:
            exons_count_all[specie] = final_fine_exon_intron_sub_dfs[(final_fine_exon_intron_sub_dfs['exon_type'] == 'fine') & (
            final_fine_exon_intron_sub_dfs['seqid'].str.startswith(specie.upper()))].shape[0]

        average_introns_length_all = round(float(filtered_introns_df['length'].mean()), 2)

        print(f'{round(((unique_seqid_in_final_df_number / total_files) * 100), 2)}% of files end up succesfully. \n '
              f' There were {len(broken_files)} errors. List of files with errors in broken_files file.'
              f'Invalid nucleotides (different than A T G C) are written with positions in invalid_nucleotides file. '
              f'\n \n Found {unique_seqid_in_final_df_number} unique files in data frame from {count_files} '
              f'from input. It means {count_files - unique_seqid_in_final_df_number} files wasn\'t written '
              f'down to output, propably because lack of exons. List of that files is in no_exons_files file. \n \n There were '
              f'{alignments_with_one_exon} alignments with only one exon in substrings, so they are not included in overall '
              f'amount of exons.')

        files_not_in_exon_intron_dfs = find_alignments_not_in_df(path_to_file_after_mafft,
                                                                 files_in_final_fine_exon_intron_sub_dfs)

        with open(f'results_{current_time}', 'w') as file:
            file.write(
                f'{round((((unique_seqid_in_final_df_number) / total_files) * 100), 2)}% of files end up succesfully. \n '
                f' There were {len(broken_files)} errors. List of files with errors in broken_files file.'
                f'Invalid nucleotides (different than A T G C) are written with positions in invalid_nucleotides file. '
                f'\n \n Found {unique_seqid_in_final_df_number} unique files in data frame from {count_files} '
                f'from input. It means {count_files - unique_seqid_in_final_df_number} files wasn\'t written '
                f'down to output, propably because lack of exons. List of that files is in no_exons_files file. \n \n There were '
                f'{alignments_with_one_exon} alignments with only one exon in substrings, so they are not included in overall '
                f'amount of exons.')

        with open('no_exons_files', 'w') as file:
            file.write('Files absent in data frame from input. Probably because lack of introns \n')
            for filename in files_not_in_exon_intron_dfs:
                file.write(filename + '\n')

    else:
        print("No valid data frames to concatenate.")
        final_fine_exon_intron_sub_dfs = pd.DataFrame(columns=column_names_substring)
        introns_count_all = 0

    # Saving table in tsv gff format
    # save_to_gff_file(all_exons_df, 'all_exons_gff.tsv')
    # save_to_gff_file(tlh_exons_df, 'tlh_exons_gff.tsv')
    # save_to_gff_file(fine_exons_df, 'fine_exons_gff.tsv')
    # save_to_gff_file(final_fine_exon_intron_sub_dfs, 'final_fine_exon_intron_sub_dfs.tsv')


    with open('broken_files', 'w') as file:
        file.write('Files which made error: \n')
        file.write(str(broken_files))

    with open('invalid_nucleotides', 'w') as file:
        for alignment_id, nt1, nt2, pos in invalid_nucleotides:
            file.write(f'{alignment_id}\t{nt1}\t{nt2}\t{pos}\n')

    return (exons_count_all, introns_count_all, unique_seqid_in_final_df_number_all, average_introns_length,
            average_introns_length_all, substrings_introns_number_all, final_fine_exon_intron_sub_dfs)




########################################################################################################################
#####################################               MINOR FUNCTIONS                #####################################
########################################################################################################################


def manual_counting_alignment_score(seq1, seq2, target_length):
    score = sum(1 for a,b in zip(seq1, seq2) if a == b)
    total = min(len(seq1), len(seq2))
    identity_level = round((score / total * 100), 2)
    return score, identity_level


def percentage_of_advancement(directory):
    count = 0
    for file in os.listdir(directory):
        if file.endswith(".fasta"):
            count += 1
    return count


def cleaning_gaps_from_both_edges(sequence_transcipt, sequence_genome, alignment, filename):
    sequence_transcipt_left_shorted = sequence_transcipt.lstrip("-")
    count_deleted_gaps_left = len(sequence_transcipt) - len(
        sequence_transcipt_left_shorted)  # inform how many gaps were deleted from 5'

    sequence_transcipt_right_and_left_shorted = sequence_transcipt_left_shorted.rstrip("-")
    count_deleted_gaps_right = len(sequence_transcipt_left_shorted) - len(
        sequence_transcipt_right_and_left_shorted)  # inform how many gaps were deleted from 3'

    if count_deleted_gaps_right:
        sequence_genome_left_and_right_shorted = sequence_genome[count_deleted_gaps_left:-count_deleted_gaps_right]
    else:
        sequence_genome_left_and_right_shorted = sequence_genome[count_deleted_gaps_left:]

    if len(sequence_genome_left_and_right_shorted) != len(sequence_transcipt_right_and_left_shorted):
        raise ValueError(
            f" VALUE ERROR: Given sequences: {alignment[0].id} - len: {len(alignment[0])} and  {alignment[1].id} -"
            f" len {len(alignment[1])} from {filename} must have the same length.")
    return sequence_transcipt_right_and_left_shorted, sequence_genome_left_and_right_shorted, count_deleted_gaps_left


def indices_of_introns_and_exons_fun(seq1_dt, seq2_dt, seq1, seq2, count_deleted_gaps_left, alignment_id,
                                     invalid_nucleotides):
    # linking pre-exon's nucleotides (nt-nt pairs) and pre-intron's nucleotides (gap-nt pairs) to theirs indices
    # NOTE: output is nucleotide or gap from seq1 (transcript sequence). Even if it is gap like this:
    # seq1: aaatttggg, seq2: aaa---ggg, output will be 'ttt' instead '---'
    # it shows that where were gap, seq2 genome or in seq1 transcript
    temp_exon = []  # list containing indices of intron's positions in sequence.
    temp_intron = []  # list containing indices of intron's positions in sequence.
    index = 0
    paired_nucleotides = zip(seq1, seq2)
    valid_nucleotides = set('ATGCatgc')

    if len(seq1_dt) != len(seq2_dt):
        print(seq1)
        print(seq2)
        raise ValueError("Two sequences must be of the same length.")

    for nt1, nt2 in paired_nucleotides:
        if "-" in (nt1, nt2):
            temp_intron.append(index + count_deleted_gaps_left)
        elif nt1 not in valid_nucleotides or nt2 not in valid_nucleotides:
            if alignment_id not in invalid_nucleotides:
                invalid_nucleotides.append(
                    (alignment_id, nt1, nt2, index + count_deleted_gaps_left))  # adding nucleotides and positions
        else:
            temp_exon.append(index + count_deleted_gaps_left)
        index += 1  # index to showing position in sequence
    return temp_intron, temp_exon, invalid_nucleotides


def start_and_end_parameters_of_exons_dict_fun(temp_exon, gaps_signs):
    # linking indices into single strand of possible exons
    start_exon = temp_exon[1]
    all_exon_range_dict = {}  # that dictionary contains: key(index of start exon) and value(index of end exon)

    for i in range(len(temp_exon)):
        # print(start_exon)
        if (temp_exon[i] - temp_exon[i - 1]) > len(gaps_signs) +1:  # if difference between two indices of exons's
            # positions is bigger than given number (2), that smaller number is index of 3' nucleotide in exon
            end_exon = temp_exon[i - 1]  # temp_exon[i-1] because it is index in list temp_exon. +1 because python
            # starts counting from 0
            all_exon_range_dict[start_exon] = end_exon + 1  # creating dictionary with all exons, even with theese too
            # short and with too low homology
            start_exon = temp_exon[i] + 1
    all_exon_range_dict[start_exon] = temp_exon[-1]  # last pair

    if all_exon_range_dict == {}:
        raise ValueError(f"ValueError - all_exon_range_dict seems to be empty. Its length = {len(all_exon_range_dict)}")
    # print(f" \n all_exon_range_dict: {all_exon_range_dict}")
    return all_exon_range_dict


def start_and_end_parameters_of_introns_dict_fun(all_exon_range_dict):
    # linking indices into single strand of possible introns
    intron_range_dict = {}
    keys_from_all_exon_range_dict = sorted(all_exon_range_dict.keys())  # list cointaining exons's start positions
    # print(f"keys_from_all_exon_range_dict: {keys_from_all_exon_range_dict}\n")

    for i in range(len(keys_from_all_exon_range_dict) - 1):
        start_intron = all_exon_range_dict[keys_from_all_exon_range_dict[i]] + 1
        end_intron = keys_from_all_exon_range_dict[i + 1] - 1
        intron_range_dict[start_intron] = end_intron
    if intron_range_dict == {}:
        raise ValueError(f"ValueError - intron_range_dict seems to be empty. Its length = {len(intron_range_dict)}")
    # print(f" \n intron_range_dict: {intron_range_dict}")
    return intron_range_dict


def making_exons_gff_table_fun(all_exon_range_dict, min_length_aligned_sequence, extreme_homology, seq1_dt, seq2_dt,
                               alignment_id, column_names, filename):
    exon_rows_to_concat_table = []

    aligner = Align.PairwiseAligner()
    aligner.mismatch_score = 0  # customized setting towards get pure percentage of homology, not alignment with gap
    # penalty etc. The object of interest is simply that how many nt in query has equivalend in target seq.
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0

    for key in all_exon_range_dict:
        query = seq1_dt[key - 1: all_exon_range_dict[key]]
        target = seq2_dt[key - 1: all_exon_range_dict[key]]
        target_length = len(target)

        if target_length >= min_length_aligned_sequence:
            score, identity_level = manual_counting_alignment_score(query, target, target_length)

            # print(f'query: {query}, \ntarget: {target}, percent of homology = {round(score / len(target), 2)* 100} %, '
            #       f'score = {score}, len = {len(target)}\n')
            if identity_level >= extreme_homology * 100:
                exon_class = 'fine'  # print(f"\n query: {query} \n target: {target} \n with start: {key} and end: {all_exon_range_dict[key]}  # goes to fine_exons with homology: {identity_level} \n")  # if len(target) != len(query):  # print(f'target: {len(target)} query: {len(query)}')
            else:
                exon_class = 'tlh'  # print(f"\n query: {query} \n target: {target} \n with start: {key} and end: {all_exon_range_dict[key]}  # goes to tlh_exons with homology: {identity_level}")  # if len(target) != len(query):  # print(f'target: {len(target)} query: {len(query)}')
        else:
            identity_level = 0
            exon_class = "ts"  # print(f"\n query: {query} \n target: {target} \n with start: {key} and end: {all_exon_range_dict[key]}  # goes to ts_exons with homology: {identity_level}")  # if len(target) != len(query):  # print(f'target: {len(target)} query: {len(query)}')
        source = extending_source_data_frame(alignment_id)
        strand = extending_strand_data_frame(alignment_id)
        exon_rows_to_concat_table.append((alignment_id, source, exon_class, key, (all_exon_range_dict[key]),
                                          ((all_exon_range_dict[key]) - key + 1), identity_level, strand, ".",
                                          filename))
    df = pd.DataFrame(exon_rows_to_concat_table, columns=column_names)

    return exon_rows_to_concat_table, df


def extending_source_data_frame(alignment_id):
    # print(f"Alignment: {alignment[0].id}, funkcja: extending_source_data_frame, przeszlo")
    keywords = ["TRINITY", "BACKBONE", "SCAFFOLD"]
    found = False
    for key in keywords:
        alignment_id = str(alignment_id).lower()
        if alignment_id.find(key.lower()) != -1:
            source = key
            found = True
            break
    if not found:
        source = None  # none means unknown
    return source


def extending_strand_data_frame(alignment_id):
    # print(f"Alignment: {alignment[0].id}, funkcja: extending_strand_data_frame, przeszlo")
    if str(alignment_id).find("SL+") != -1:
        strand = "+"
    elif str(alignment_id).find("SL-") != -1:
        strand = "-"
    else:
        strand = None
    return strand


def adding_introns_to_gff_data_frame(all_exons_df, tlh_exons_df, fine_exons_df, column_names):
    data_frames_list = [all_exons_df, tlh_exons_df, fine_exons_df]

    for i, df in enumerate(data_frames_list):
        df.sort_values(by=['seqid', 'start'], inplace=True)

        df['seqid_introns'] = df['seqid'].shift(
            -1)  # creating new column with seqid of next sequention. In next steps, rows without matching seqid and
        # seqid_introns, would not be consideresd
        df['intron_start'] = df['start'].shift(-1)  # creating index of start position in introns

        intron_df_temp = df[df['seqid'] == df['seqid_introns']].copy()
        intron_df_temp['start'] = df['end'] + 1
        intron_df_temp['end'] = df['intron_start'] - 1
        intron_df_temp['exon_type'] = 'intron'
        intron_df_temp['homology'] = 0
        intron_df_temp['length'] = intron_df_temp['end'] - intron_df_temp['start'] + 1
        intron_df_temp = intron_df_temp[column_names]

        df = df[column_names]
        df = pd.concat([df, intron_df_temp]).sort_values(by=['seqid', 'start']).reset_index(drop=True)

        data_frames_list[i] = df

    return data_frames_list[0], data_frames_list[1], data_frames_list[2]


def adding_introns_to_final_gff_data_frame(df, column_names_substring):
    df.sort_values(by=['seqid', 'start'], inplace=True)
    df['seqid_introns'] = df['seqid'].shift(-1)  # creating new column with seqid of next sequention.
    # In next steps, rows without matching seqid and seqid_introns, would not be consideresd
    df['intron_start'] = df['start'].shift(-1)  # creating index of start position in introns

    intron_df_temp = df[df['seqid'] == df['seqid_introns']].copy()
    intron_df_temp['start'] = df['end'] + 1
    intron_df_temp['end'] = df['intron_start'] - 1
    intron_df_temp['exon_type'] = 'intron'
    intron_df_temp['homology'] = 0
    intron_df_temp['length'] = intron_df_temp['end'] - intron_df_temp['start'] + 1
    intron_df_temp = intron_df_temp[column_names_substring]

    df = df[column_names_substring]
    df = pd.concat([df, intron_df_temp]).sort_values(by=['seqid', 'start']).reset_index(drop=True)

    return df


def save_to_gff_file(gff_final_data_frame, filename):
    # print(f"Alignment: {alignment[0].id}, funkcja: save_to_gff_file, przeszlo")
    if os.path.isfile(filename):
        user_input = input(f"GFF file  already exists. Do you want overwrite? y/n: ")
        if user_input == "n":
            base_name, ext = os.path.splitext(filename)
            i = 1
            while os.path.isfile(f"{base_name}_{i}{ext}"):
                i += 1
            filename = f"{base_name}_{i}{ext}"
            gff_final_data_frame.to_csv(filename, sep="\t", index=False)

        elif user_input == "y":
            gff_final_data_frame.to_csv(filename, sep="\t", index=False)
        else:
            print("Clarify your answer. Nothing has done.")

    else:
        gff_final_data_frame.to_csv(filename, sep="\t", index=False)


def extracting_strands_from_alignment(
        alignment):  # that function describe what strand is transcript strand and genomic strand - in our files,
    # transcript strand has longer ID.
    # print(f"Alignment: {alignment[0].id}, funkcja: extracting_strands_from_alignment, przeszlo")
    seq1 = alignment[0].seq
    seq2 = alignment[1].seq
    # print(f"\n Just run function: extracting_strands_from_alignment")
    return seq1, seq2  # return seq1_is_transcript #obecnie nie uzywane


def getting_files_in_directory(path_to_file_after_mafft):
    with os.scandir(path_to_file_after_mafft) as entries:
        return [entry.name.split('.')[0] for entry in entries if entry.is_file()]


def getting_files_from_df(df, column_name, files_in_final_fine_exon_intron_sub_dfs_all, species):
    # creating list with unique names in df and counting unique names in df
    unique_names = df[column_name].astype(str).unique().tolist()
    for specie in species:
        files_in_final_fine_exon_intron_sub_dfs_all[specie] = [file for file in unique_names if
                                                               file.upper().startswith(str(specie).upper())]
    unique_numbers = df[column_name].nunique()
    return files_in_final_fine_exon_intron_sub_dfs_all, unique_names, unique_numbers


def find_alignments_not_in_df(path_to_file_after_mafft, files_in_exons_df):
    files_in_directory = set(getting_files_in_directory(path_to_file_after_mafft))
    files_in_exons_df = set(files_in_exons_df)
    files_not_in_df = files_in_directory - files_in_exons_df
    return files_not_in_df


def check_exon_homology(start, end, seq1_dt, seq2_dt, min_length_aligned_sequence, extreme_homology):
    query = seq1_dt[start - 1: end]
    target = seq2_dt[start - 1: end]
    target_length = len(target)
    # print(f'przed if: \n query: target: \n{query}\n{target} \n target length: {target_length}\n minlength: {min_length_aligned_sequence}')
    if target_length >= min_length_aligned_sequence:
        score, identity_level = manual_counting_alignment_score(query, target, target_length)
        if identity_level >= extreme_homology * 100:
            # print(f'score: {score}, identity level: {identity_level}, extreme homology: {extreme_homology}')
            return True, score


def create_substring_dict(all_exon_range_dict, seq1_dt, seq2_dt, min_length_aligned_sequence, extreme_homology,
                          alignments_with_one_exon):
    fine_exon_sub_dict = {}
    fine_intron_sub_dict = {}
    current_exon_substring = []
    current_intron_substring = {}
    substring_index = 1

    for start, end in all_exon_range_dict.items():
        if check_exon_homology(start, end, seq1_dt, seq2_dt, min_length_aligned_sequence, extreme_homology):
            current_exon_substring.append((start, end))

        else:
            if len(current_exon_substring) >= 2:  #it means if there are two exons what equals there is intron between
                fine_exon_sub_dict[f"substring_{substring_index}"] = dict(current_exon_substring)
                substring_index += 1
                current_exon_substring = []
            else:
                alignments_with_one_exon += 1
                current_exon_substring = []

    if len(current_exon_substring) >= 2:  # to add last substring when last nucleotides in alignment are exon.
        fine_exon_sub_dict[f"substring_{substring_index}"] = dict(current_exon_substring)

    fine_intron_sub_dict = create_intron_substring_dict(fine_exon_sub_dict, fine_intron_sub_dict)

    return fine_exon_sub_dict, fine_intron_sub_dict, alignments_with_one_exon


def create_intron_substring_dict(fine_exon_sub_dict, fine_intron_sub_dict):
    for substring, exon_ranges in fine_exon_sub_dict.items():
        current_intron_substring = {}
        sorted_start_exon_ranges = sorted(exon_ranges.keys())
        for i in range(len(sorted_start_exon_ranges) - 1):
            start = exon_ranges[sorted_start_exon_ranges[i]] + 1
            end = sorted_start_exon_ranges[i + 1] - 1
            current_intron_substring[start] = end
        fine_intron_sub_dict[substring] = current_intron_substring
    return fine_intron_sub_dict


def create_final_exon_gff_table_from_substring(fine_exon_sub_dict, seq1_dt, seq2_dt, alignment_id,
                                               column_names_substring, filename):
    fine_exon_sub_df = []
    exon_class = 'fine'
    source = extending_strand_data_frame(alignment_id)
    strand = extending_strand_data_frame(alignment_id)
    for substring, exons_range in fine_exon_sub_dict.items():
        for start, end in exons_range.items():
            query = seq1_dt[start - 1: end]
            target = seq2_dt[start - 1: end]
            target_length = len(target)
            score, identity_level = manual_counting_alignment_score(query, target, target_length)
            GC = GC_counting(start, end, seq1_dt)
            fine_exon_sub_df.append((
                alignment_id, source, exon_class, start, end, (end - start + 1), identity_level, GC, strand, '.', filename,
                substring))
    data_frame = pd.DataFrame(fine_exon_sub_df, columns=column_names_substring)
    return data_frame


def create_final_intron_table_from_substring(fine_intron_sub_dict, alignment_id, column_names_substring, filename, seq2_dt):
    fine_intron_sub_df = []
    exon_class = 'intron'
    source = extending_strand_data_frame(alignment_id)
    strand = extending_strand_data_frame(alignment_id)

    for substring, exon_range in fine_intron_sub_dict.items():
        for start, end in exon_range.items():
            GC = GC_counting(start, end, seq2_dt)
            fine_intron_sub_df.append(
                (alignment_id, source, exon_class, start, end, (end - start +1), 0, GC, strand, '.', filename, substring))
    data_frame = pd.DataFrame(fine_intron_sub_df, columns=column_names_substring)
    return data_frame

def GC_counting(start, end, seq):
    sequence = seq[start - 1: end]
    # print(sequence)
    GC_count = round((sequence.count('g') + sequence.count('c')) / len(sequence) * 100, 2)
    # print(GC_count)
    return GC_count

def extracting_informations_to_return(df_with_introns, list_of_species, introns_count_all, average_introns_length,
                                      substrings_introns_number_all, unique_seqid_in_final_df_number_all):
    for specie in list_of_species:
        specie_df = df_with_introns[df_with_introns['seqid'].str.contains(str(specie).upper())]
        introns_count_all[specie] = round(len(specie_df), 2)
        average_introns_length[specie] = round(float(specie_df['length'].mean()), 2)
        substrings_introns_number_all[specie] = (specie_df[['seqid', 'substring']].drop_duplicates()).shape[0]
        unique_seqid_in_final_df_number_all[specie] = specie_df['seqid'].nunique()

    return introns_count_all, average_introns_length, substrings_introns_number_all, unique_seqid_in_final_df_number_all


if __name__ == '__main__':
    cutting_scrap(path_to_file_before_mafft, path_to_file_after_mafft, 2, 0.93, 20)

