#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import subprocess
import os
import sys
# sys.path.append('/home/norbert/mrdn/euglena/kod_i_pliki/Nowa_wersja_funkcji')
sys.path.append('/home/szala/euglena/kod_i_pliki')
from nowa_wersja_funkcji import cutting_scrap



path_to_file_before_mafft = '/home/szala/euglena/kod_i_pliki/divided_fastas'
path_to_file_after_mafft = '/home/szala/euglena/kod_i_pliki/fastas_after_mafft'

# path_to_file_before_mafft = ('/home/norbert/mrdn/euglena/kod_i_pliki/surowe_pliki_plus_minus_500/raw_reads_9/'
#                              'merging_fastas')
# path_to_file_after_mafft = ('/home/norbert/mrdn/euglena/kod_i_pliki/surowe_pliki_plus_minus_500/raw_reads_9/'
#                             'fastas_after_mafft_na_probe')

# path_to_save_charts = '/home/norbert/mrdn/euglena/kod_i_pliki/Nowa_wersja_funkcji/charts_homology'
path_to_save_charts = '/home/szala/euglena/kod_i_pliki/charts_homology'



def ten_times_funkcja():
    exons_count_dict = {}
    fine_introns_count_dict = {}
    unique_seqid_in_final_df_number_all_dict = {}
    average_introns_length_dict = {}
    substrings_introns_number_dict = {}
    average_introns_length_dict_all = {}

    for i in range(90, 98, 1):
        i = int(i)/100
        (exons_count_all, fine_introns_count, unique_seqid_in_final_df_number_all, average_introns_length,
         average_introns_length_all, substrings_introns_number, main_df) = (cutting_scrap(path_to_file_before_mafft,
                                                                                          path_to_file_after_mafft, 2,
                                                                                          i, 20))

        exons_count_dict[i] = exons_count_all
        exons_count_df = making_summary_df(exons_count_dict)
        exons_count_df['description'] = 'exons_count'

        fine_introns_count_dict[i] = fine_introns_count
        fine_introns_count_df = making_summary_df(fine_introns_count_dict)
        fine_introns_count_df['description'] = 'fine_introns_count'

        unique_seqid_in_final_df_number_all_dict[i] = unique_seqid_in_final_df_number_all
        unique_seqid_in_final_df_number_all_df = making_summary_df(unique_seqid_in_final_df_number_all_dict)
        unique_seqid_in_final_df_number_all_df['description'] = 'unique_seqid_in_final_df_number_all'

        average_introns_length_dict[i] = average_introns_length
        average_introns_length_dict_all[i] = average_introns_length_all

        average_introns_length_df_spec = pd.DataFrame(average_introns_length_dict).T.reset_index().rename(
            columns={'index': 'min_ex_homology'})
        average_introns_length_df_all = pd.DataFrame(list(average_introns_length_dict_all.items()),
                                                     columns=['min_ex_homology', 'average_introns_length'])
        average_introns_length_df = pd.merge(average_introns_length_df_spec, average_introns_length_df_all,
                                             on='min_ex_homology', how='outer')
        average_introns_length_df['description'] = 'average_introns_length'
        average_introns_length_df = average_introns_length_df.rename(columns={'average_introns_length':'total'})

        substrings_introns_number_dict[i] = substrings_introns_number
        substrings_introns_number_df = making_summary_df(substrings_introns_number_dict)
        substrings_introns_number_df['description'] = 'substrings_introns_number'

    enormous_df = pd.concat([fine_introns_count_df, unique_seqid_in_final_df_number_all_df, average_introns_length_df,
                            substrings_introns_number_df], ignore_index=True)

    return enormous_df, exons_count_df, fine_introns_count_df, unique_seqid_in_final_df_number_all_df, average_introns_length_df, substrings_introns_number_df


def making_summary_df(dictionary):
    print(dictionary)
    df = pd.DataFrame(dictionary).T.reset_index().rename(columns={'index': 'min_ex_homology'})
    df['sum'] = df[['gra', 'lon', 'hie']].sum(axis=1)
    return df


if __name__ == '__main__':
    (enormous_df, exons_count_df, fine_introns_count_df, unique_seqid_in_final_df_number_all_df,
     average_introns_length_df, substrings_introns_number_df) = ten_times_funkcja()


# Wykresy

#zmienne ogólne:
df_count = fine_introns_count_df
df_unique = unique_seqid_in_final_df_number_all_df
df_sub = substrings_introns_number_df
df_exons = exons_count_df
df_av_int_len = average_introns_length_df

bar_width = 0.2
opacity = 0.8


# Liczba intronow vs minimalna długość
fig, ax = plt.subplots(figsize=(12, 8))
index = np.arange(len(df_count))

# Zmiana kolejności słupków i ich pozycji
bars_sum = ax.bar(index - 1.5*bar_width, df_count['sum'], bar_width, alpha=opacity, edgecolor='black', color='#0d1949',
                  label='Sum')
bar_gra = ax.bar(index - 0.5*bar_width, df_count['gra'], bar_width, alpha=opacity, edgecolor='black', color='#3e5a6f',
                 label='E. gracilis')
bar_lon = ax.bar(index + 1.5*bar_width, df_count['lon'], bar_width, alpha=opacity, edgecolor='black', color='#a0ddba',
                 label='E. longa')
bar_hie = ax.bar(index + 0.5*bar_width, df_count['hie'], bar_width, alpha=opacity, edgecolor='black', color='#6f9c94',
                 label='E. hiemalis')

ax.grid(True, linestyle='--', alpha=0.7)

# Dodanie wartości nad słupkami
for bars in [bars_sum, bar_gra, bar_lon, bar_hie]:
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval, round(yval, 2), ha='center', va='bottom', fontsize=8)

ax.set_ylabel('Number of introns', fontsize=12, labelpad=15)
ax.set_xlabel('Minimal exon\'s homology', fontsize=12, labelpad=15)
ax.set_title('Number of introns depending on minimal exon\'s homology', fontsize=18, pad=20)

ax.set_xticks(index)
ax.set_xticklabels(df_count['min_ex_homology']*100)
ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_facecolor('#f5f5f5')

handles = [bars_sum, bar_gra, bar_hie, bar_lon]
ax.legend(handles=handles, labels=['Sum', 'E. gracilis', 'E. hiemalis', 'E. longa'], loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4)

plt.tight_layout()
plt.savefig(f'{path_to_save_charts}/number_of_introns.png', format='png', dpi=300)
plt.show()




###
# Liczba intronów na sekwencję vs minimalna długość
fig, ax = plt.subplots(figsize=(12, 8))
index = np.arange(len(df_unique))

# Zmiana kolejności słupków i ich pozycji
bars_sum = ax.bar(index - 1.5*bar_width, df_count['sum']/df_unique['sum'], bar_width, alpha=opacity, edgecolor='black', color='#0d1949',
                  label='Sum')
bar_gra = ax.bar(index - 0.5*bar_width, df_count['gra']/df_unique['gra'], bar_width, alpha=opacity, edgecolor='black', color='#3e5a6f',
                 label='E. gracilis')
bar_lon = ax.bar(index + 1.5*bar_width, df_count['lon']/df_unique['lon'], bar_width, alpha=opacity, edgecolor='black', color='#a0ddba',
                 label='E. longa')
bar_hie = ax.bar(index + 0.5*bar_width, df_count['hie']/df_unique['hie'], bar_width, alpha=opacity, edgecolor='black', color='#6f9c94',
                 label='E. hiemalis')

ax.grid(True, linestyle='--', alpha=0.7)

# Dodanie wartości nad słupkami
for bars in [bars_sum, bar_gra, bar_lon, bar_hie]:
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval, f"{yval:.2f}", ha='center', va='bottom', fontsize=8)
ax.set_ylabel('Number of introns per sequence', fontsize=12, labelpad=15)
ax.set_xlabel('Minimal exon\'s homology', fontsize=12, labelpad=15)
ax.set_title('Number of introns per sequence depending on minimal exon\'s homology', fontsize=18, pad=20)

ax.set_xticks(index)
ax.set_xticklabels(df_unique['min_ex_homology']*100)
ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_facecolor('#f5f5f5')

ax.legend(handles=handles, labels=['Sum', 'E. gracilis', 'E. hiemalis', 'E. longa'], loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4)

plt.tight_layout()
plt.savefig(f'{path_to_save_charts}/number_of_introns_per_sequence.png', format='png', dpi=300)
plt.show()
###



###
# Liczba intronów na substring vs minimalna długość
fig, ax = plt.subplots(figsize=(12, 8))
index = np.arange(len(df_unique))

# Zmiana kolejności słupków i ich pozycji
bars_sum = ax.bar(index - 1.5*bar_width, df_count['sum']/df_sub['sum'], bar_width, alpha=opacity, edgecolor='black', color='#0d1949',
                  label='Sum')
bar_gra = ax.bar(index - 0.5*bar_width, df_count['gra']/df_sub['gra'], bar_width, alpha=opacity, edgecolor='black', color='#3e5a6f',
                 label='E. gracilis')
bar_lon = ax.bar(index + 1.5*bar_width, df_count['lon']/df_sub['lon'], bar_width, alpha=opacity, edgecolor='black', color='#a0ddba',
                 label='E. longa')
bar_hie = ax.bar(index + 0.5*bar_width, df_count['hie']/df_sub['hie'], bar_width, alpha=opacity, edgecolor='black', color='#6f9c94',
                 label='E. hiemalis')

ax.grid(True, linestyle='--', alpha=0.7)

# Dodanie wartości nad słupkami
for bars in [bars_sum, bar_gra, bar_lon, bar_hie]:
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval, f"{yval:.2f}", ha='center', va='bottom', fontsize=8)
ax.set_ylabel('Number of introns per substring', fontsize=12, labelpad=15)
ax.set_xlabel('Minimal exon\'s homology', fontsize=12, labelpad=15)
ax.set_title('Number of introns per substring depending on minimal exon\'s homology', fontsize=18, pad=20)

ax.set_xticks(index)
ax.set_xticklabels(df_sub['min_ex_homology']*100)
ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_facecolor('#f5f5f5')

ax.legend(handles=handles, labels=['Sum', 'E. gracilis', 'E. hiemalis', 'E. longa'], loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4)


plt.tight_layout()
plt.savefig(f'{path_to_save_charts}/number_of_introns_per_substring.png', format='png', dpi=300)
plt.show()
###


###
# Liczba intronów na ekson vs minimalna długość - kwestionuje zasadność tego porównania ale dlaczego nie?
# (teoreteycznie intronów powinno być tyle co n eksonów - 1, ale tutaj zmieniamy warunki gry bo odrzucamy niektóre
# eksony  i introny więc ten stosunek powinien być mniejszy
fig, ax = plt.subplots(figsize=(12, 8))
index = np.arange(len(df_exons))

# Zmiana kolejności słupków i ich pozycji
bars_sum = ax.bar(index - 1.5*bar_width, df_exons['sum']/df_count['sum'], bar_width, alpha=opacity, edgecolor='black', color='#0d1949',
                  label='Sum')
bar_gra = ax.bar(index - 0.5*bar_width, df_exons['gra']/df_count['gra'], bar_width, alpha=opacity, edgecolor='black', color='#3e5a6f',
                 label='E. gracilis')
bar_lon = ax.bar(index + 1.5*bar_width, df_exons['lon']/df_count['lon'], bar_width, alpha=opacity, edgecolor='black', color='#a0ddba',
                 label='E. longa')
bar_hie = ax.bar(index + 0.5*bar_width, df_exons['hie']/df_count['hie'], bar_width, alpha=opacity, edgecolor='black', color='#6f9c94',
                 label='E. hiemalis')

ax.grid(True, linestyle='--', alpha=0.7)

# Dodanie wartości nad słupkami
for bars in [bars_sum, bar_gra, bar_lon, bar_hie]:
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval, f"{yval:.2f}", ha='center', va='bottom', fontsize=8)
ax.set_ylabel('Number of exons per intron', fontsize=12, labelpad=15)
ax.set_xlabel('Minimal exon\'s homology', fontsize=12, labelpad=15)
ax.set_title('Number of exons per intron depending on minimal exon\'s homology', fontsize=18, pad=20)

ax.set_xticks(index)
ax.set_xticklabels(df_exons['min_ex_homology']*100)
ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_facecolor('#f5f5f5')

ax.legend(handles=handles, labels=['Sum', 'E. gracilis', 'E. hiemalis', 'E. longa'], loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4)

plt.tight_layout()
plt.savefig(f'{path_to_save_charts}/number_of_exons_per_intron.png', format='png', dpi=300)
plt.show()
###




###
# Średnia długość intronów vs minimalna długość
fig, ax = plt.subplots(figsize=(12, 8))
index = np.arange(len(df_av_int_len))

# Zmiana kolejności słupków i ich pozycji
bars_sum = ax.bar(index - 1.5*bar_width, df_av_int_len['total'], bar_width, alpha=opacity, edgecolor='black', color='#0d1949',
                  label='Total')
bar_gra = ax.bar(index - 0.5*bar_width, df_av_int_len['gra'], bar_width, alpha=opacity, edgecolor='black', color='#3e5a6f',
                 label='E. gracilis')
bar_lon = ax.bar(index + 1.5*bar_width, df_av_int_len['lon'], bar_width, alpha=opacity, edgecolor='black', color='#a0ddba',
                 label='E. longa')
bar_hie = ax.bar(index + 0.5*bar_width, df_av_int_len['hie'], bar_width, alpha=opacity, edgecolor='black', color='#6f9c94',
                 label='E. hiemalis')

ax.grid(True, linestyle='--', alpha=0.7)

# Dodanie wartości nad słupkami
for bars in [bars_sum, bar_gra, bar_lon, bar_hie]:
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval, f"{yval:.2f}", ha='center', va='bottom', fontsize=8)
ax.set_ylabel('Average intron\'s length', fontsize=12, labelpad=15)
ax.set_xlabel('Minimal exon\'s homology', fontsize=12, labelpad=15)
ax.set_title('Average intron\'s length depending on minimal exon\'s homology', fontsize=18, pad=20)

ax.set_xticks(index)
ax.set_xticklabels(df_av_int_len['min_ex_homology']*100)
ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_facecolor('#f5f5f5')

ax.legend(handles=handles, labels=['Total', 'E. gracilis', 'E. hiemalis', 'E. longa'], loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4)

plt.tight_layout()
plt.savefig(f'{path_to_save_charts}/average_introns_length.png', format='png', dpi=300)
plt.show()
###



# liczba unikalnych sekwencji vs minimalna długość
fig, ax = plt.subplots(figsize=(12, 8))
index = np.arange(len(df_unique))

# Zmiana kolejności słupków i ich pozycji
bars_sum = ax.bar(index - 1.5*bar_width, df_unique['sum'], bar_width, alpha=opacity, edgecolor='black', color='#0d1949',
                  label='Sum')
bar_gra = ax.bar(index - 0.5*bar_width, df_unique['gra'], bar_width, alpha=opacity, edgecolor='black', color='#3e5a6f',
                 label='E. gracilis')
bar_lon = ax.bar(index + 1.5*bar_width, df_unique['lon'], bar_width, alpha=opacity, edgecolor='black', color='#a0ddba',
                 label='E. longa')
bar_hie = ax.bar(index + 0.5*bar_width, df_unique['hie'], bar_width, alpha=opacity, edgecolor='black', color='#6f9c94',
                 label='E. hiemalis')

ax.grid(True, linestyle='--', alpha=0.7)

# Dodanie wartości nad słupkami
for bars in [bars_sum, bar_gra, bar_lon, bar_hie]:
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval, f"{yval:.2f}", ha='center', va='bottom', fontsize=8)
ax.set_ylabel('Number of unique ID sequences', fontsize=12, labelpad=15)
ax.set_xlabel('Minimal exon\'s homology', fontsize=12, labelpad=15)
ax.set_title('Number of unique ID sequences depending on minimal exon\'s homology', fontsize=18, pad=20)

ax.set_xticks(index)
ax.set_xticklabels(df_unique['min_ex_homology']*100)
ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_facecolor('#f5f5f5')

ax.legend(handles=handles, labels=['Sum', 'E. gracilis', 'E. hiemalis', 'E. longa'], loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4)

plt.tight_layout()
plt.savefig(f'{path_to_save_charts}/number_of_unique_seqid.png', format='png', dpi=300)
plt.show()
###




###
# liczba substringów vs minimalna długość
fig, ax = plt.subplots(figsize=(12, 8))
index = np.arange(len(df_sub))

# Zmiana kolejności słupków i ich pozycji
bars_sum = ax.bar(index - 1.5*bar_width, df_sub['sum'], bar_width, alpha=opacity, edgecolor='black', color='#0d1949',
                  label='Sum')
bar_gra = ax.bar(index - 0.5*bar_width, df_sub['gra'], bar_width, alpha=opacity, edgecolor='black', color='#3e5a6f',
                 label='E. gracilis')
bar_lon = ax.bar(index + 1.5*bar_width, df_sub['lon'], bar_width, alpha=opacity, edgecolor='black', color='#a0ddba',
                 label='E. longa')
bar_hie = ax.bar(index + 0.5*bar_width, df_sub['hie'], bar_width, alpha=opacity, edgecolor='black', color='#6f9c94',
                 label='E. hiemalis')

ax.grid(True, linestyle='--', alpha=0.7)

# Dodanie wartości nad słupkami
for bars in [bars_sum, bar_gra, bar_lon, bar_hie]:
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval, f"{yval:.2f}", ha='center', va='bottom', fontsize=8)
ax.set_ylabel('Number of substrings ', fontsize=12, labelpad=15)
ax.set_xlabel('Minimal exon\'s homology', fontsize=12, labelpad=15)
ax.set_title('Number of substrings depending on minimal exon\'s homology', fontsize=18, pad=20)

ax.set_xticks(index)
ax.set_xticklabels(df_sub['min_ex_homology']*100)
ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_facecolor('#f5f5f5')

ax.legend(handles=handles, labels=['Sum', 'E. gracilis', 'E. hiemalis', 'E. longa'], loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4)
plt.tight_layout()
plt.savefig(f'{path_to_save_charts}/number_of_substring.png', format='png', dpi=300)
plt.show()
###


###
# liczba substringów na sekwencje vs minimalna długość
fig, ax = plt.subplots(figsize=(12, 8))
index = np.arange(len(df_sub))

# Zmiana kolejności słupków i ich pozycji
bars_sum = ax.bar(index - 1.5*bar_width, df_sub['sum']/df_unique['sum'], bar_width, alpha=opacity, edgecolor='black', color='#0d1949',
                  label='Sum')
bar_gra = ax.bar(index - 0.5*bar_width, df_sub['gra']/df_unique['gra'], bar_width, alpha=opacity, edgecolor='black', color='#3e5a6f',
                 label='E. gracilis')
bar_lon = ax.bar(index + 1.5*bar_width, df_sub['lon']/df_unique['lon'], bar_width, alpha=opacity, edgecolor='black', color='#a0ddba',
                 label='E. longa')
bar_hie = ax.bar(index + 0.5*bar_width, df_sub['hie']/df_unique['hie'], bar_width, alpha=opacity, edgecolor='black', color='#6f9c94',
                 label='E. hiemalis')

ax.grid(True, linestyle='--', alpha=0.7)

# Dodanie wartości nad słupkami
for bars in [bars_sum, bar_gra, bar_lon, bar_hie]:
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval, f"{yval:.2f}", ha='center', va='bottom', fontsize=8)
ax.set_ylabel('Number of substrings per sequence', fontsize=12, labelpad=15)
ax.set_xlabel('Minimal exon\'s homology', fontsize=12, labelpad=15)
ax.set_title('Number of substrings per sequence depending on minimal exon\'s homology', fontsize=18, pad=20)

ax.set_xticks(index)
ax.set_xticklabels(df_sub['min_ex_homology']*100)
ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_facecolor('#f5f5f5')

ax.legend(handles=handles, labels=['Sum', 'E. gracilis', 'E. hiemalis', 'E. longa'], loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4)

plt.tight_layout()
plt.savefig(f'{path_to_save_charts}/number_of_substring_per_sequence.png', format='png', dpi=300)
plt.show()
###




########################################### EXONY ####################################
###
# Liczba exonow vs minimalna długość
fig, ax = plt.subplots(figsize=(12, 8))
index = np.arange(len(df_exons))

# Zmiana kolejności słupków i ich pozycji
bars_sum = ax.bar(index - 1.5*bar_width, df_exons['sum'], bar_width, alpha=opacity, edgecolor='black', color='#0d1949',
                  label='Sum')
bar_gra = ax.bar(index - 0.5*bar_width, df_exons['gra'], bar_width, alpha=opacity, edgecolor='black', color='#3e5a6f',
                 label='E. gracilis')
bar_lon = ax.bar(index + 1.5*bar_width, df_exons['lon'], bar_width, alpha=opacity, edgecolor='black', color='#a0ddba',
                 label='E. longa')
bar_hie = ax.bar(index + 0.5*bar_width, df_exons['hie'], bar_width, alpha=opacity, edgecolor='black', color='#6f9c94',
                 label='E. hiemalis')

ax.grid(True, linestyle='--', alpha=0.7)

# Dodanie wartości nad słupkami
for bars in [bars_sum, bar_gra, bar_lon, bar_hie]:
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval, f"{yval:.2f}", ha='center', va='bottom', fontsize=8)
ax.set_ylabel('Number of exons', fontsize=12, labelpad=15)
ax.set_xlabel('Minimal exon\'s homology', fontsize=12, labelpad=15)
ax.set_title('Number of exons depending on minimal exon\'s homology', fontsize=18, pad=20)

ax.set_xticks(index)
ax.set_xticklabels(df_exons['min_ex_homology']*100)
ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_facecolor('#f5f5f5')

ax.legend(handles=handles, labels=['Sum', 'E. gracilis', 'E. hiemalis', 'E. longa'], loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4)

plt.tight_layout()
plt.savefig(f'{path_to_save_charts}/number_of_exons.png', format='png', dpi=300)
plt.show()
###



# Liczba exonow na sekwencję vs minimalna długość
fig, ax = plt.subplots(figsize=(12, 8))
index = np.arange(len(df_exons))

# Zmiana kolejności słupków i ich pozycji
bars_sum = ax.bar(index - 1.5*bar_width, df_exons['sum']/df_unique['sum'], bar_width, alpha=opacity, edgecolor='black', color='#0d1949',
                  label='Sum')
bar_gra = ax.bar(index - 0.5*bar_width, df_exons['gra']/df_unique['gra'], bar_width, alpha=opacity, edgecolor='black', color='#3e5a6f',
                 label='E. gracilis')
bar_lon = ax.bar(index + 1.5*bar_width, df_exons['lon']/df_unique['lon'], bar_width, alpha=opacity, edgecolor='black', color='#a0ddba',
                 label='E. longa')
bar_hie = ax.bar(index + 0.5*bar_width, df_exons['hie']/df_unique['hie'], bar_width, alpha=opacity, edgecolor='black', color='#6f9c94',
                 label='E. hiemalis')

ax.grid(True, linestyle='--', alpha=0.7)

# Dodanie wartości nad słupkami
for bars in [bars_sum, bar_gra, bar_lon, bar_hie]:
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval, f"{yval:.2f}", ha='center', va='bottom', fontsize=8)
ax.set_ylabel('Number of exons per sequence', fontsize=12, labelpad=15)
ax.set_xlabel('Minimal exon\'s homology', fontsize=12, labelpad=15)
ax.set_title('Number of exons per sequence depending on minimal exon\'s homology', fontsize=18, pad=20)

ax.set_xticks(index)
ax.set_xticklabels(df_unique['min_ex_homology']*100)
ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_facecolor('#f5f5f5')

ax.legend(handles=handles, labels=['Sum', 'E. gracilis', 'E. hiemalis', 'E. longa'], loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4)

plt.tight_layout()
plt.savefig(f'{path_to_save_charts}/number_of_exons_per_sequence.png', format='png', dpi=300)
plt.show()
###



###
# Liczba exonow na substring vs minimalna długość
fig, ax = plt.subplots(figsize=(12, 8))
index = np.arange(len(df_exons))

# Zmiana kolejności słupków i ich pozycji
bars_sum = ax.bar(index - 1.5*bar_width, df_exons['sum']/df_sub['sum'], bar_width, alpha=opacity, edgecolor='black', color='#0d1949',
                  label='Sum')
bar_gra = ax.bar(index - 0.5*bar_width, df_exons['gra']/df_sub['gra'], bar_width, alpha=opacity, edgecolor='black', color='#3e5a6f',
                 label='E. gracilis')
bar_lon = ax.bar(index + 1.5*bar_width, df_exons['lon']/df_sub['lon'], bar_width, alpha=opacity, edgecolor='black', color='#a0ddba',
                 label='E. longa')
bar_hie = ax.bar(index + 0.5*bar_width, df_exons['hie']/df_sub['hie'], bar_width, alpha=opacity, edgecolor='black', color='#6f9c94',
                 label='E. hiemalis')

ax.grid(True, linestyle='--', alpha=0.7)

# Dodanie wartości nad słupkami
for bars in [bars_sum, bar_gra, bar_lon, bar_hie]:
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval, f"{yval:.2f}", ha='center', va='bottom', fontsize=8)
ax.set_ylabel('Number of exons per substring', fontsize=12, labelpad=15)
ax.set_xlabel('Minimal exon\'s homology', fontsize=12, labelpad=15)
ax.set_title('Number of exons per substring depending on minimal exon\'s homology', fontsize=18, pad=20)

ax.set_xticks(index)
ax.set_xticklabels(df_sub['min_ex_homology']*100)
ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_facecolor('#f5f5f5')

ax.legend(handles=handles, labels=['Sum', 'E. gracilis', 'E. hiemalis', 'E. longa'], loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4)

plt.tight_layout()
plt.savefig(f'{path_to_save_charts}/number_of_exons_per_substring.png', format='png', dpi=300)
plt.show()
###