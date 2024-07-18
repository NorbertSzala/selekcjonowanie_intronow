#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
# sys.path.append('/home/norbert/mrdn/euglena/kod_i_pliki/Nowa_wersja_funkcji')
sys.path.append('/home/szala/euglena/kod_i_pliki')
from nowa_wersja_funkcji import cutting_scrap


# na serwer
path_to_file_before_mafft = "/home/szala/euglena/kod_i_pliki/divided_fastas"
path_to_file_after_mafft = "/home/szala/euglena/kod_i_pliki/fastas_after_mafft"


# lokalnie
# path_to_file_before_mafft = ('/home/norbert/mrdn/euglena/kod_i_pliki/surowe_pliki_plus_minus_500/raw_reads_9/'
#                              'merging_fastas')
# path_to_file_after_mafft = ('/home/norbert/mrdn/euglena/kod_i_pliki/surowe_pliki_plus_minus_500/raw_reads_9/'
#                             'fastas_after_mafft_na_probe')

def comparing_df():
    dfs = {}
    for i in range(90, 92, 1):
        results = (cutting_scrap(path_to_file_before_mafft, path_to_file_after_mafft, 2, int(i)/100, 20))
        print(len(results))
        (exons_count_all, introns_count_all, unique_seqid_in_final_df_number_all, average_introns_length,
         average_introns_length_all, substrings_introns_number_all) = results[:6]

        dfs[f'df_{i}'] = results[6]
    return dfs


dfs = comparing_df()
df90 = dfs['df_90']
df90 = df90.iloc[:, :-1]
df91 = dfs['df_91']
df91 = df91.iloc[:, :-1]

df_diff = df90.merge(df91, how='left', indicator=True).query('_merge == "left_only"').drop(columns=['_merge'])

added_exons = (df_diff['exon_type'] == 'intron').sum()
print(added_exons)

df_diff.to_csv('df_diff_homology.tsv', sep="\t", index=False)

