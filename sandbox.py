from pyfaidx import Fasta
import multiprocessing as mp
import pandas as pd
import time
from collections import defaultdict
import kmertools as kt
import os

def prepare_directory(parent='../results/'):
    if parent is None:
        parent = './results/'
    import datetime
    directory = parent + datetime.datetime.now().strftime("results_%d%b%Y-%H%M/")
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory


def zip_chrom(key, df, ref):
    fa = Fasta(ref)
    seq = str(fa[key])
    varray = [[] for _ in range(len(seq))]
    for index, row in df.iterrows():
        varray[row['POS'] - 1].append((row['ALT'], row['REF']))
    return [key, zip(list(seq), varray)]


def write_ref_dict_file(ref_dict):
    directory = prepare_directory(parent='./ref_var_dict/')
    for key, value in ref_dict:
        fp = open(directory + key + '.csv', 'w')
        fp.write('REF_fasta\tALT\tREF_vcf')
        for rec in value:
            for i in rec:
                fp.write(str(i) + '\t')
        fp.close()


def prepare_reference_dict(fasta_path, variants, delim='\t'):
    """
    :param primary_chroms: boolean True means only include original autosomal chromosomes, False includes everything
    :param delim: indicates how your variant file is separated
    :param fasta_path: path to file containing reference sequence
    :param variants: path to bed file containing the variants in the reference sequence
    :return: a dictionary mapping chromosome names to an array of tuples containing the reference allele in the first index and the variant allele in the second index if it exists
    """
    start = time.time()
    var_df = pd.read_csv(variants, sep=delim)
    keys = kt.get_primary_chroms(fasta_path)[7:8]
    args = []
    for key in keys:
        args.append((key, var_df[var_df.iloc[:, 0] == key], fasta_path))
    pool = mp.Pool(mp.cpu_count())
    results = [funccall.get() for funccall in [pool.starmap_async(zip_chrom, args)]]
    pool.close()
    print('Done processing variants in %f' % (time.time() - start))
    directory = prepare_directory(parent='./ref_var_dict/')
    for chrom_key in results[0]:
        fp = open(directory + chrom_key[0] + '.csv', 'w')
        fp.write('REF_fasta\tALT\tREF_vcf')
        for rec in chrom_key[1]:
            for i in rec:
                fp.write(str(i) + '\t')
            fp.write('\n')
        fp.close()
    print('Reference dictionary prepared in %f' % (time.time() - start))
    return


if __name__ == "__main__":
    f, v, vb = kt.test_data()
    prepare_reference_dict(f, vb)
    print('Done')
