import pandas as pd
import matplotlib.pyplot as plt
import kmertools as kt

transitions_path = '/Users/simonelongo/Documents/QuinlanLabFiles/kmer_data/results/kp_21oct2019/bp_counts_per3mer.csv'
counts_path = '/Users/simonelongo/Documents/QuinlanLabFiles/kmer_data/data/ref_genome_kmer_freq.csv'

counts = pd.read_csv(counts_path, index_col=0)
counts.columns = ['count']
transitions = pd.read_csv(transitions_path, index_col=0)
merged = transitions.join(counts, how='inner')

freq = merged.iloc[:, 0:4].div(merged['count'], axis=0)

flank_count_AG = pd.DataFrame()
flank_count_CT = pd.DataFrame()
for idx, row in freq.iterrows():
    if idx[1] == 'A':
        flank_count_AG.loc[idx[0], idx[2]] = row['G']
    if idx[1] == 'T':
        flank_count_AG.loc[idx[0], idx[2]] = row['C']
    if idx[1] == 'G':
        flank_count_CT.loc[idx[0], idx[2]] = row['A']
    if idx[1] == 'C':
        flank_count_CT.loc[idx[0], idx[2]] = row['T']


def plot_flanks(flank_count, title='K-mer transitions'):
    flank_count = flank_count.sort_index(ascending=False)
    flank_count = flank_count.reindex(sorted(flank_count.columns), axis=1)
    plt.matshow(flank_count)
    plt.colorbar()
    plt.title(title)
    plt.xticks(range(4), list(flank_count.columns))
    plt.yticks(range(4), list(flank_count.index))


plot_flanks(flank_count_AG)
plot_flanks(flank_count_CT)
"""
G > T == C > A
G > A == C > T
G > C == C > G
T > A == A > T
T > C == A > G  # one from paper I want to reproduce
T > G == A > C
"""
