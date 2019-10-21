from cyvcf2 import VCF
from pyfaidx import Fasta
from collections import defaultdict, Counter
import multiprocessing as mp
import pandas as pd
import kmertools as kt
import time

REF_FASTA = ""
VCF_PATH = ""
VARIANT_CSV = ""
KMER_SIZE = 3


def process_region(region):
    vcf = VCF(VCF_PATH)
    ref = Fasta(REF_FASTA)
    transitions = defaultdict(Counter)
    variant_positions = defaultdict(kt.Variant)
    start_idx_offset = int(KMER_SIZE / 2 + 1)
    kmer_mid_idx = int(start_idx_offset - 1)  # also halfway index for kmer
    for section in region:
        print('Processing ' + str(section))
        start = time.time()
        for variant in vcf(str(section)):
            if kt.is_quality_variant(variant):
                new_var = kt.Variant(variant.REF, "".join(variant.ALT), variant.POS, variant.CHROM)
                variant_positions[variant.POS] = new_var
                # take 7mer around variant. pyfaidx excludes start index and includes end index
                adj_seq = ref[str(variant.CHROM)][(variant.POS - start_idx_offset):(variant.POS + kmer_mid_idx)].seq
                if kt.complete_sequence(adj_seq):
                    transitions[adj_seq][variant.ALT[0]] += 1
        print('Time to complete section ' + str(section) + ': ' + str(time.time() - start))
    return [transitions, variant_positions]


def run_vcf_analysis(nprocs):
    regions = kt.get_split_vcf_regions(VCF_PATH, nprocs)
    pool = mp.Pool(nprocs)
    results = [funccall.get() for funccall in [pool.map_async(process_region, regions)]]
    pool.close()
    return results


def vcf_singleton_analysis(vcf_path, ref_fasta, kmer_length, nprocs=0, var_out_file='genome_variants.bed', store=True):
    if nprocs == 0:
        nprocs = mp.cpu_count()
    global REF_FASTA, VCF_PATH, VARIANT_CSV, KMER_SIZE
    REF_FASTA = ref_fasta
    KMER_SIZE = kmer_length
    if kt.is_vcf(vcf_path):
        transitions_list = []
        VCF_PATH = vcf_path
        results = run_vcf_analysis(nprocs)

        output = open(var_out_file, "a+")
        output.write("CHROM\tPOS\tREF\tALT\n")  # header same for all
        for result in results[0]:
            for collection in result:
                if isinstance(list(collection.values())[0], kt.Variant):
                    for k, v in collection.items():
                        output.write(str(v))
                else:
                    transitions_list.append(collection)
        output.close()
        merged_dict = kt.merge_defaultdict(transitions_list)
        if store:
            count_fpath = "bp_counts_per" + str(KMER_SIZE) + "mer.csv"
            transitions_df = pd.DataFrame.from_dict(merged_dict, orient='index')
            transitions_df.to_csv(count_fpath)
        return merged_dict
    else:  # VCF has already been processed, process text file to save time
        VARIANT_CSV = vcf_path
