import os
import random
from collections import Counter, defaultdict
import pandas as pd
from pyfaidx import Fasta
from cyvcf2 import VCF, Writer
import itertools



def append_variants_to_vcf(chrom, start, stop):
    return "tabix /Users/simonelongo/too_big_for_icloud/gnomad.genomes.r2.1.1.sites.vcf.bgz " + str(chrom) + ":" + str(
        start) + "-" + str(
        stop) + " >> samp.vcf"


def generate_sample_vcf(filename='/Users/simonelongo/too_big_for_icloud/gnomad.genomes.r2.1.1.sites.vcf.bgz'):
    """Takes a large VCF file and takes random samples from each chromosome to make a smaller VCF for testing"""
    reference_genome = Fasta('/Users/simonelongo/too_big_for_icloud/REFERENCE_GENOME_GRch37.fa')
    vcf = VCF(filename)
    write = Writer('samp.vcf', vcf)
    write.write_header()
    chrom_num = list(reference_genome.keys())
    for chrom in chrom_num[0:25]:
        if chrom != 'MT':
            begin = random.randint(1000, len(reference_genome[chrom]) - 1000)
            os.system(append_variants_to_vcf(chrom, begin, begin + 1000))
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


def find_ref_kmer_freq(kmer_length):
    """Count of occurrences of each kmer on the reference genome"""
    expected_path = 'data/chr22_' + str(kmer_length) + 'mer_frequency.csv'
    col_names = ['ref_count']
    if os.path.exists(expected_path):
        df = pd.read_csv(expected_path, header=None, index_col=0, names=col_names)
        print('Reference kmer frequency successfully read.')
        return df
    counts = Counter()
    ref_genome = Fasta('/Users/simonelongo/too_big_for_icloud/REFERENCE_GENOME_GRch37.fa')
    ref_seq = str(ref_genome["22"])
    for i in range(len(ref_seq) - (kmer_length - 1)):  # This takes the (1-based) reference sequence for chromosome 22
        next_seq = ref_seq[i:(i + kmer_length)]
        if not ('N' in next_seq or 'n' in next_seq):
            counts[next_seq] += 1
    print('Reference kmer frequency population done.')
    outfile = pd.DataFrame.from_dict(counts, orient='index')
    outfile.to_csv(expected_path)
    outfile = pd.read_csv(expected_path, header=None, skiprows=1, index_col=0, names=col_names)

    return outfile


def is_vcf(f_path):
    """Returns True if filepath is a valid VCF file"""
    tokens = os.path.splitext(f_path)
    if tokens[1] == '.vcf' or tokens[1] == '.gz':
        return True
    return False


def is_quality_variant(var_to_test):
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
        if is_quality_variant(variant):
            # join is required because 'ALT' is returned as a list
            variant_positions[variant.CHROM][variant.POS] = Variant(variant.REF, "".join(variant.ALT), variant.POS, variant.CHROM)
    return variant_positions
