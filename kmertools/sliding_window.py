import kmertools as kt
import os
import multiprocessing as mp
import pandas as pd


def window_search(f_path, kmersize):
    df = pd.read_csv(f_path, sep=r'\s+')
    observed_vars = 0
    for index, row in df.iterrows():
        if row['ALT'] is not None:
            observed_vars += 1

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

