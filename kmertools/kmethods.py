import os
import sys
import random
import time
from collections import Counter, defaultdict
import pandas as pd
import multiprocessing as mp
import numpy as np
from pyfaidx import Fasta
from cyvcf2 import VCF, Writer
import itertools
from kmertools import Variant, Storage, Kmer


def ran_seq(chunk):
    sequence = ""
    for i in range(chunk):
        sequence += random.choice('ACTG')
    return sequence


def gen_random_sequence(length):
    pool = mp.Pool()
    nprocs = mp.cpu_count()
    chunk_size = int(length / nprocs) + 1
    args = []
    tot = chunk_size
    while tot <= length:
        args.append(chunk_size)
        tot += chunk_size
    new_tot = np.array(args).sum()
    if new_tot < length:
        args.append(length - new_tot)
    print(args)

    results = [funccall.get() for funccall in [pool.map_async(ran_seq, args)]]
    pool.close()
    random_seq = ""
    for seq in results[0]:
        random_seq += seq
    return random_seq


def append_variants_to_vcf(chrom, start, stop):
    return "tabix /Users/simonelongo/too_big_for_icloud/gnomad.genomes.r2.1.1.sites.vcf.bgz " + str(chrom) + ":" + str(
        start) + "-" + str(
        stop) + " >> samp.vcf"


# TODO: Fix so that reference is not used
def generate_sample_vcf(filename='/Users/simonelongo/too_big_for_icloud/gnomad.genomes.r2.1.1.sites.vcf.bgz'):
    """Takes a large VCF file and takes random samples from each chromosome to make a smaller VCF for testing"""
    vcf = VCF(filename)
    write = Writer('samp.vcf', vcf)
    write.write_header()
    for chrom_num, chrom_len in get_primary_chroms(filename):
        begin = random.randint(1000, chrom_len - 1000)
        os.system(append_variants_to_vcf(chrom_num, begin, begin + 1000))
    write.close()


def complement(c):
    base_pairs = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
        'N': 'N'
    }
    try:
        return base_pairs[c.upper()]
    except KeyError:
        raise ValueError(c + " is not a valid nucleotide.")


def get_complementary_sequence(sequence):
    """
    Returns a string of nucleotides complementary to the input string
    All letters in input sequence must A, C, T, or G, otherwise will raise a ValueError
    """
    comp_seq = []
    for c in sequence[::-1]:  # take the reverse complement
        comp_seq.append(complement(c))
    return "".join(comp_seq)


def generate_kmers(k):
    """Generates a list of all possible DNA sequences of length k. E.g. generate_kmers(2) will return
    [ AA, AC, AT, AG, CA, CC, CT, CG, TA, TC, TT, TG, GA, GC, GT, GG ] """
    len_k = int(k)
    if len_k < 0:
        raise ValueError("Must be a positive integer")
    combos = list(itertools.product('ACTG', repeat=len_k))
    seqs = []
    for seq in combos:
        seqs.append(''.join(seq))
    return seqs


def find_ref_kmer_freq(kmer_length, ref_fasta):
    """Count of occurrences of each kmer on the reference genome"""
    expected_path = 'data/chr22_' + str(kmer_length) + 'mer_frequency.csv'
    col_names = ['ref_count']
    if os.path.exists(expected_path):
        df = pd.read_csv(expected_path, header=None, index_col=0, names=col_names)
        print('Reference kmer frequency successfully read.')
        return df
    counts = Counter()
    ref_genome = Fasta(ref_fasta)
    ref_seq = ""
    for chrom in ref_genome.keys():
        ref_seq += str(ref_genome[chrom])
    for i in range(len(ref_seq) - (kmer_length - 1)):  # This takes the (1-based) reference sequence for chromosome 22
        next_seq = ref_seq[i:(i + kmer_length)]
        if not ('N' in next_seq or 'n' in next_seq):
            counts[next_seq] += 1
    print('Reference kmer frequency population done.')
    outfile = pd.DataFrame.from_dict(counts, orient='index')
    outfile.to_csv(expected_path)
    outfile = pd.read_csv(expected_path, header=None, skiprows=1, index_col=0, names=col_names)

    return outfile


def ref_genome_as_string(ref_fasta, keys=None):
    ref_genome = Fasta(ref_fasta)
    if keys is None:
        keys = ref_genome.keys()
    ref_seq = ""
    for chrom in keys:
        ref_seq += str(ref_genome[chrom])
    return ref_seq


def is_vcf(f_path):
    """Returns True if filepath is a valid VCF file"""
    tokens = os.path.splitext(f_path)
    allowed_extensions = ['.vcf', '.gz', '.bgz']
    return tokens[1] in allowed_extensions


def is_quality_variant(var_to_test):
    """
    high quality variants will have FILTER == None
    AND we are ignoring insertions and deltions here
    """
    return var_to_test.FILTER is None and len(var_to_test.ALT) == 1 \
           and len(var_to_test.REF) == 1 and len(var_to_test.ALT[0]) == 1


def is_quality_nonsingleton(var_to_test):
    """
    high quality variants will have FILTER == None
    Additionally, variants shoud NOT be singletons ('AC' != 1) meaning that it is a unique observation
    AND we are ignoring insertions and deltions here
    """
    return var_to_test.FILTER is None and var_to_test.INFO.get('AC') != 1 and len(var_to_test.ALT) == 1 \
           and len(var_to_test.REF) == 1 and len(var_to_test.ALT[0]) == 1


def is_quality_singleton(var_to_test):
    """
    high quality variants will have FILTER == None
    Additionally, variants shoud be singletons ('AC' == 1) meaning that it is a unique observation
    AND we are ignoring insertions and deltions here
    """
    return var_to_test.FILTER is None and var_to_test.INFO.get('AC') == 1 and len(var_to_test.ALT) == 1 \
           and len(var_to_test.REF) == 1 and len(var_to_test.ALT[0]) == 1


def complete_sequence(adj_seq):
    """If 'N' is present iin the sequence, kmer is undefined"""
    return not ('N' in adj_seq or 'n' in adj_seq)


def get_transition(ref, alt):
    """Gives a common scheme for transitions when considering reverse complements"""
    if ref != 'A' and ref != 'C':
        ref = complement(ref)
        alt = complement(alt)
    return ref + alt


def get_variants(filepath):
    """
    Returns a defaultdict containing all singleton variants contained in the input VCF file:
    Default dict uses chromosome number as a key and maps to a defaultdict of variants arranged by position
    """
    from kmertools.kclass import Variant
    variant_positions = defaultdict(lambda: defaultdict(Variant))
    for variant in VCF(filepath):
        if is_quality_singleton(variant):
            # join is required because 'ALT' is returned as a list
            variant_positions[variant.CHROM][variant.POS] = Variant(variant.REF, "".join(variant.ALT), variant.POS,
                                                                    variant.CHROM)
    return variant_positions


def process_variants_old(variants, kmer_size, ref_fasta):
    print('This method has been deprecated... for now.')
    if not kmer_size % 2 == 1:
        kmer_size += 1  # kmer size must be an odd number.
    if not isinstance(variants, pd.DataFrame):
        print("Reading variant file.")
        variants_df = pd.read_csv(variants, sep='\t', low_memory=False)
    ref = Fasta(ref_fasta)
    transitions = defaultdict(Counter)
    start_idx_offset = int(kmer_size / 2 + 1)
    kmer_mid_idx = int(start_idx_offset - 1)  # also halfway index for kmer
    print("Processing variant file")
    width = 40
    setup_progressbar(toolbar_width=width)
    df_prog_inc = int(len(variants_df) / width)
    for idx, row in variants_df.iterrows():
        if idx % df_prog_inc == 0:
            sys.stdout.write("-")
            sys.stdout.flush()
        position = row['POS']
        # take 7mer around variant. pyfaidx excludes start index and includes end index
        adj_seq = ref[str(row['CHROM'])][(position - start_idx_offset):(position + kmer_mid_idx)].seq
        if complete_sequence(adj_seq):
            transitions[adj_seq][row['ALT']] += 1
    sys.stdout.write("]\n")  # this ends the progress bar

    count_fpath = "bp_counts_per" + str(kmer_size) + "mer.csv"
    transitions_df = pd.DataFrame.from_dict(transitions, orient='index')
    # merged_df = pd.merge(transitions_df, find_ref_kmer_freq(kmer_size), left_index=True, right_index=True, how='inner')
    transitions_df.to_csv(count_fpath)
    return transitions_df


def generate_csv_from_variants(variants, outfile="variants_samp.csv", close=False):
    """Converts a dictionary to a csv to prevemt redundant slow operations"""
    # if not type(variants) == dict:
    #     print("Input must be a dictionary")
    #     return
    output = open(outfile, "a+")
    if not os.path.exists(outfile):
        output.write("POS\tREF\tALT\n")  # header same for all
    for k, v in variants.items():
        output.write(str(v))
    if close:
        output.close()


def split_pandas_rows(df):
    """Taken from https://stackoverflow.com/questions/40357434/pandas-df-iterrows-parallelization"""
    # create as many processes as there are CPUs on your machine
    num_processes = mp.cpu_count()

    # calculate the chunk size as an integer
    chunk_size = int(df.shape[0] / num_processes)

    # this solution was reworked from the above link.
    # will work even if the length of the dataframe is not evenly divisible by num_processes
    chunks = [df.ix[df.index[i:i + chunk_size]] for i in range(0, df.shape[0], chunk_size)]


def generate_heatmap(transitions, kmer_length):
    kmer_freq = find_ref_kmer_freq(kmer_length)
    merged = transitions.join(kmer_freq, how='inner')
    # TODO: implement


def setup_progressbar(toolbar_width=40):
    """taken from: https://stackoverflow.com/questions/3160699/python-progress-bar"""
    # setup toolbar
    sys.stdout.write("[%s]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width + 1))  # return to start of line, after '['


def process_vcf_region(vcf_path, regions):
    vcf = VCF(vcf_path)
    variant_positions = defaultdict(Variant)
    for region in regions:
        for variant in vcf(str(region)):
            if is_quality_singleton(variant):
                variant_positions[variant.POS] = Variant(variant.REF, "".join(variant.ALT), variant.POS, variant.CHROM)
    print("VCF chunk processed")
    return variant_positions


# FIXME
# def process_vcf_region2(regions, kmer_size=3):
#     """
#     Process VCF Regions 2: Electric Boogaloo
#         Returns a list of variants and the transitions in one sweep. That's wassup!
#     """
#     vcf = VCF(constants.VCF_PATH)
#     # variant_positions = defaultdict(Variant)
#     ref = Fasta(REF_FASTA_PATH)
#     transitions = defaultdict(Counter)
#     variant_positions = defaultdict(Variant)
#     start_idx_offset = int(kmer_size / 2 + 1)
#     kmer_mid_idx = int(start_idx_offset - 1)  # also halfway index for kmer
#     for region in regions:
#         for variant in vcf(str(region)):
#             if is_quality_variant(variant):
#                 new_var = Variant(variant.REF, "".join(variant.ALT), variant.POS, variant.CHROM)
#                 variant_positions[variant.POS] = new_var
#                 # take 7mer around variant. pyfaidx excludes start index and includes end index
#                 adj_seq = ref[str(variant.CHROM)][(variant.POS - start_idx_offset):(variant.POS + kmer_mid_idx)].seq
#                 if complete_sequence(adj_seq):
#                     transitions[adj_seq][variant.ALT[0]] += 1
#     print("VCF chunk processed")
#     return [transitions, variant_positions]


# def run_vcf_parallel(vcf_path, nprocs):
#     if constants.VCF_PATH != vcf_path:
#         constants.VCF_PATH = vcf_path
#     regions = get_split_vcf_regions(vcf_path, nprocs)
#     pool = mp.Pool(nprocs)
#     results = [funccall.get() for funccall in [pool.map_async(process_vcf_region2, regions)]]
#     pool.close()
#     return results


# def vcf_parallel(vcf_path, nprocs, outfile='genome_variants_agg.csv', store=True):
#     if nprocs == 0:
#         nprocs = mp.cpu_count()
#     results = run_vcf_parallel(vcf_path, nprocs)
#     output = open(outfile, "a+")
#     output.write("CHROM\tPOS\tREF\tALT\n")  # header same for all
#     transitions_list = []
#     for result in results[0]:
#         if isinstance(list(result[0].values())[0], Variant):
#             for k, v in result.items():
#                 output.write(str(v))
#         else:
#             transitions_list.append(result[0])
#     output.close()
#     print(outfile + " saved.")
#     print("Done reading VCF file.")
#     merged_dict = merge_defaultdict(transitions_list)
#     if store:
#         count_fpath = "bp_counts_per" + str(constants.KMER_SIZE) + "mer.csv"
#         transitions_df = pd.DataFrame.from_dict(merged_dict, orient='index')
#         transitions_df.to_csv(count_fpath)
#     return merged_dict


def get_split_vcf_regions(vcf_path, nprocs):
    vcf = VCF(vcf_path)
    from kmertools import VCFRegion
    num_entries = np.sum(vcf.seqlens)
    chunk_size = int(num_entries / nprocs) + 1
    names = get_primary_chroms(vcf_path)
    num_chunk = 0
    regions = []
    gen_pos = 0
    current_chromosome = 0
    chrom_pos = 0
    while num_chunk < nprocs and current_chromosome < len(names):
        current_chunk = 0
        region = []
        while current_chunk < chunk_size:
            if current_chromosome >= len(names):
                current_chunk = chunk_size
                continue
            remaining_chunk = chunk_size - current_chunk
            if remaining_chunk <= (vcf.seqlens[current_chromosome] - chrom_pos):
                new_region = VCFRegion(vcf.seqnames[current_chromosome], chrom_pos, chrom_pos + remaining_chunk)
                region.append(new_region)
                chrom_pos += remaining_chunk
                current_chunk += new_region.size()
                gen_pos += new_region.size()
                continue
            else:  # remaining chunk can fit remainder of chromosome and then some
                new_region = VCFRegion(vcf.seqnames[current_chromosome], chrom_pos, vcf.seqlens[current_chromosome])
                region.append(new_region)
                current_chunk += new_region.size()
                gen_pos += new_region.size()
                chrom_pos = 0
                current_chromosome += 1
        regions.append(region)
    # print(regions)
    return regions


def merge_defaultdict(dict_list):
    master_count = defaultdict(Counter)
    for counts in dict_list:
        for k, v in counts.items():
            for alt, count in v.items():
                master_count[k][alt] += count
    return master_count


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


def kmer_search(sequence, kmer_length):
    counts = Counter()
    for i in range(len(sequence) - (kmer_length - 1)):  # This takes the (1-based) reference sequence for chromosome 22
        next_seq = sequence[i:(i + kmer_length)]
        if not ('N' in next_seq or 'n' in next_seq):
            counts[next_seq] += 1
    return counts


def get_kmer_count(sequence, kmer_length):
    start = time.time()
    args = split_seq(sequence, mp.cpu_count(), overlap=kmer_length)
    args = [[seq, kmer_length] for seq in args]
    pool = mp.Pool(mp.cpu_count())
    results = [res.get() for res in [pool.starmap_async(kmer_search, args)]]
    pool.close()
    counts = Counter()
    for result in results[0]:
        for k, v in result.items():
            counts[k] += v
    print("Done in " + str(time.time() - start))
    return counts


def test_data(big=False):
    fasta = '/Users/simonelongo/too_big_for_icloud/ref_genome/REFERENCE_GENOME_GRch37.fa'
    big_vcf = '/Users/simonelongo/too_big_for_icloud/gnomad.genomes.r2.1.1.sites.vcf.bgz'
    lil_vcf = '/Users/simonelongo/Documents/QuinlanLabFiles/kmer_data/data/samp.vcf.bgz'
    if big:
        return fasta, big_vcf
    else:
        return fasta, lil_vcf


def combine_complements(dataframe):
    # Pandas data frame indexed by k-mer sequence, with numeric values
    total = defaultdict(Counter)
    for index, row in dataframe.iterrows():
        seq = Kmer(index)
        for i in row.index:
            # print(str(seq), str(row[i]))
            total[seq][i] += row[i]
    return pd.DataFrame.from_dict(total, orient='index')


def get_primary_chroms(f_path):
    split = os.path.splitext(f_path)
    tokens = []
    while len(split[1]) > 0:
        tokens.append(split[1])
        split = os.path.splitext(split[0])
    if '.fa' in tokens:
        fa = Fasta(f_path)
        chrom_keys = [n for n in fa.keys() if len(n) < 6 and 'chrM' not in n]
    if '.vcf' in tokens:
        v = VCF(f_path)
        chrom_keys = [n for n in v.seqnames if len(n) < 6 and 'chrM' not in n]
    return chrom_keys


def prepare_directory(parent='../results/'):
    if parent is None:
        parent = '../results/'
    import datetime
    directory = parent + datetime.datetime.now().strftime("results_%d%b%Y-%H%M/")
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory
