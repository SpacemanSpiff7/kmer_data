import kmertools as kt
import os
import multiprocessing as mp
import pandas as pd
import numpy as np
from collections import Counter


def count_kmers(kmer_length, ref_seq):
    counts = Counter()
    for i in range(len(ref_seq) - (kmer_length - 1)):  # This takes the (1-based) reference sequence for chromosome 22
        next_seq = ref_seq[i:(i + kmer_length)]
        if not ('N' in next_seq or 'n' in next_seq):
            counts[next_seq] += 1
    return counts


def calculate_expected(seq, alts, kmersize):
    if kmersize == 3:
        ref_count = pd.read_csv('/Users/simonelongo/Documents/QuinlanLabFiles/kmer_data/data/ref_kmer_freqs/norm_hg382_3mer_frequency.csv').to_dict()
        pass
    elif kmersize == 5:
        pass
    elif kmersize == 7:
        pass
    elif kmersize == 9:
        pass
    else:
        print('Please provide reference %d-mer count'%kmersize)
        return -1, -1
    kmers = []
    midpt = int(kmersize/2)
    observed = 0
    for start in range(len(seq) - kmersize + 1):
        next_k = seq[start:start + kmersize]
        if alts[start + midpt] in ['A','C','T','G']:
            observed += 1

    return 0, 0


def split_seq(sequence, nprocs, overlap=None):
    chunk_size = int(len(sequence) / nprocs) + 1
    args = []
    start = 0
    end = chunk_size
    for proc in range(nprocs):
        if overlap is not None:
            args.append(sequence[start:(end + overlap - 1)])
        else:
            args.append(sequence[start:end])
        start = end
        end += chunk_size
        if end > len(sequence):
            end = len(sequence)
    return args


def window_search(f_path, kmersize, window_size):
    df = pd.read_csv(f_path, sep=r'\s+')
    ref = df.iloc[:, 0].to_numpy()
    refs = "".join(df.iloc[:, 0]).upper()
    alt = df.iloc[:, 1].to_numpy()
    alts = "".join(df.iloc[:, 1]).upper()
    observed_vars = 0
    expected = 0
    start = 0
    end = window_size
    while end <= len(refs):
        next_seq = refs[start:end]
        nind = next_seq.rfind('N')
        if nind == -1:
            ex, obs = calculate_expected(next_seq, alts[start:end], kmersize)
            start += 1
            end = start + window_size
        else:
            start = nind + 1
            end = start + window_size
    pass


class WindowSearch:
    def __init__(self, *args, **kwargs):
        ref_directory = ""
        if 'kmer' not in kwargs.keys():
            print('K-mer size not specified. Please enter \'kmer=<SOME INTEGER>\' in WindowSearch constructor')
            return
        self.kmer_size = kwargs['kmer']
        if 'nprocs' in kwargs.keys():
            self.nprocs = kwargs['nprocs']
        else:
            self.nprocs = mp.cpu_count()
        if all(e in ['fasta', 'vcf'] for e in kwargs.keys()):
            ref_directory = kt.prepare_reference_dict(kwargs['fasta'], kwargs['vcf'])
        if 'directory' in kwargs.keys():
            ref_directory = kwargs['directory']
            csvs = []
        for root, dirs, files in os.walk(ref_directory):
            for file in files:
                if file.endswith(".csv"):
                    csvs.append(os.path.join(root, file))
        pool = mp.Pool(self.nprocs)
        results = pool.starmap_async(window_search, zip(args, [self.kmer_size for _ in csvs]))


if __name__ == "__main__":
    directory = '/Users/simonelongo/Documents/QuinlanLabFiles/kmer_data/ref_var_dict/results_30Oct2019-1154'
    ws = WindowSearch(kmer=3, nprocs=6, directory=directory)
