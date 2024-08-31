#!/usr/bin/env python3

import pandas as pd

# file_path = '/home/norbert/mrdn/euglena/kod_i_pliki/Nowa_wersja_funkcji/final_fine_exon_intron_sub_dfs.tsv'
file_path = '/home/szala/euglena/kod_i_pliki/final_fine_exon_intron_sub_dfs.tsv'

original_df = pd.read_csv(file_path, sep='\t')


organisms = ['GRA', 'LON', 'HIE']
intron_results = {}
exon_results = {}

for organism in organisms:
    df_intron = original_df[(original_df['attributes'].str.startswith(organism)) & (original_df['exon_type'] == 'intron')]
    df_exon = original_df[(original_df['attributes'].str.startswith(organism)) & (original_df['exon_type'] == 'fine')]

    min_len_ex = int(df_exon['length'].min())
    min_len_in = int(df_intron['length'].min())

    max_len_ex = int(df_exon['length'].max())
    max_len_in = int(df_intron['length'].max())

    GC_mean_ex = float(round(df_exon['GC'].mean(), 2))
    GC_mean_in = float(round(df_intron['GC'].mean(), 2))

    intron_results[organism] = min_len_in, max_len_in, GC_mean_in
    exon_results[organism] = min_len_ex, max_len_ex, GC_mean_ex
    print(intron_results[organism])
    print(exon_results[organism])

intron_df = pd.DataFrame.from_dict(intron_results)
exon_df = pd.DataFrame.from_dict(exon_results)
print(intron_df)
print(exon_df)