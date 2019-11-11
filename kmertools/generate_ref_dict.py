from pyfaidx import Fasta
import numpy as numpy
import pandas as pd
from cyvcf2 import VCF
import time
import os
import multiprocessing as mp
import sys
import gc


def get_primary_chroms_grch38(f_path, autosomes=True):
    # v = VCF(f_path)
    if autosomes:
        return ['chr1',
                'chr2',
                'chr3',
                'chr4',
                'chr5',
                'chr6',
                'chr7',
                'chr8',
                'chr9',
                'chr10',
                'chr11',
                'chr12',
                'chr13',
                'chr14',
                'chr15',
                'chr16',
                'chr17',
                'chr18',
                'chr19',
                'chr20',
                'chr21',
                'chr22']


chpc_directory = '/uufs/chpc.utah.edu/common/home/u0319040/longo_scratch/output/lp_shared_6nov2019/'


def prepare_directory(parent=chpc_directory, new_folder=None):
    if parent is None:
        parent = '../results/'
    import datetime
    if new_folder is not None:
        if new_folder[-1] not in '/':
            new_folder += '/'
        directory = parent + new_folder
    else:
        directory = parent + datetime.datetime.now().strftime("results_%d%b%Y-%H%M/")
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory


def write_ref_var_file(directory, key, seqvar):
    fp = open(directory + key + '.tsv', 'w')
    fp.write('REF_fasta\tALT\tREF_vcf\n')
    for rec in seqvar:
        for i in rec:
            fp.write(str(i) + '\t')
        fp.write('\n')
    fp.close()
    return


def zip_chrom(key, df, ref, directory):
    start = time.time()
    print('Processing %s' % key, flush=True)
    fa = Fasta(ref)
    seq = str(fa[key])
    varray = ['n' for _ in range(len(seq))]
    for index, row in df.iterrows():
        varray[row['POS'] - 1] = ("".join(row['ALT']), row['REF'])
    write_ref_var_file(directory, key, zip(list(seq), varray))
    seq.clear()
    varray.clear()
    gc.collect()
    print('Processed %s in %f.' % (key, (time.time() - start)), flush=True)
    return  # [key, zip(list(seq), varray)]


def process_chrom(key, bed_dir, ref, directory):
    start = time.time()
    fp = open(directory + key + '.tsv', 'w')
    fp.write('REF_fasta\tALT\tREF_vcf\n')
    print('Processing %s' % key, flush=True)
    fa = Fasta(ref)
    seq = str(fa[key])
    current_pos = 0
    df = pd.read_csv((bed_dir + key + '.bed'), sep='\t')
    for index, row in df.iterrows():
        var_pos = row['POS'] - 1
        for i in range((var_pos - current_pos - 1)):
            fp.write(seq[i] + '\t' + '0' + '\t' + '0' + '\n')
        fp.write(seq[var_pos] + '\t' + "".join(row['ALT']) + '\t' + row['REF'] + '\n')
        current_pos = var_pos
    fp.close()
    print('Processed %s in %f.' % (key, (time.time() - start)), flush=True)
    return


def prepare_reference_dict(fasta_path, variants, delim='\t', primary_chroms=True, nprocs=6):
    """
    :param nprocs: number of CPUs to use
    :param primary_chroms: boolean True means only include original autosomal chromosomes, False includes everything
    :param delim: indicates how your variant file is separated
    :param fasta_path: path to file containing reference sequence
    :param variants: path to bed file containing the variants in the reference sequence
    :return: prints files to a specified directory, 1 per chromosome
    """
    start = time.time()
    fa = Fasta(fasta_path)
    final_dir = '/uufs/chpc.utah.edu/common/home/u0319040/longo_scratch/output/chroms/'
    # var_df = pd.read_csv(variants, sep=delim, low_memory=False)
    if primary_chroms:
        keys = get_primary_chroms_grch38(fasta_path)
    else:
        keys = fa.keys()
    args = []
    directory = prepare_directory(new_folder='./ref_var_dict/')
    for key in keys:
        args.append((key, final_dir, fasta_path, directory))
    pool = mp.Pool(nprocs)
    results = [funccall.get() for funccall in [pool.starmap_async(process_chrom, args)]]
    pool.close()
    print('Done processing variants in %f' % (time.time() - start), flush=True)
    return directory


if __name__ == "__main__":
    n_procs, vcf, fasta = None, None, None
    for a in sys.argv:
        try:
            n_procs = int(a)
        except ValueError:
            if 'bed' in a:
                vcf = a
            if 'fa' in a:
                fasta = a
    if None in [n_procs, vcf, fasta]:
        print('Missing Arguments\n\tNumber of CPUs:\t%s\n\tVCF path:\t%s\n\tFASTA path:\t%s\n' % (
            str(n_procs), str(vcf), str(fasta)), flush=True)
        exit(0)

    prepare_reference_dict(fasta, vcf, nprocs=n_procs)
