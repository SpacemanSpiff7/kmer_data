from cyvcf2 import VCF
from pyfaidx import Fasta
import kmertools as kt
import multiprocessing as mp
from collections import defaultdict, Counter
import time
import pandas as pd


def process_region(region, vcf_path, fasta_path, kmer_size):
    vcf = VCF(vcf_path)
    ref = Fasta(fasta_path)
    singleton_transitions = defaultdict(Counter)
    notsingleton_transitions = defaultdict(Counter)
    all_transitions = defaultdict(Counter)
    singletons = defaultdict(kt.Variant)
    not_singletons = defaultdict(kt.Variant)
    all_vars = defaultdict(kt.Variant)
    start_idx_offset = int(kmer_size / 2 + 1)
    kmer_mid_idx = int(start_idx_offset - 1)  # also halfway index for kmer
    for section in region:
        print('Processing ' + str(section))
        start_time = time.time()
        for variant in vcf(str(section)):
            if kt.is_quality_variant(variant):
                new_var = kt.Variant(variant.REF, "".join(variant.ALT), variant.POS, variant.CHROM)
                all_vars[variant.POS] = new_var
                seq_context = ref[str(variant.CHROM)][
                              (variant.POS - start_idx_offset):(variant.POS + kmer_mid_idx)].seq
                is_complete_sequence = False
                if kt.complete_sequence(seq_context):
                    is_complete_sequence = True
                    all_transitions[seq_context][variant.ALT[0]] += 1
                if kt.is_quality_singleton(variant):
                    singletons[variant.POS] = new_var
                    if is_complete_sequence:
                        singleton_transitions[seq_context][variant.ALT[0]] += 1
                if kt.is_quality_nonsingleton(variant):
                    not_singletons[variant.POS] = new_var
                    if is_complete_sequence:
                        for alt in variant.ALT:
                            notsingleton_transitions[seq_context][alt] += 1
        print('Time to complete section ' + str(section) + ': ' + str(time.time() - start_time))
        return {'all': [all_transitions, all_vars],
                'singletons': [singleton_transitions, singletons],
                'not_singletons': [notsingleton_transitions, not_singletons]}


def merge_and_save(results, destination):
    """input is a list of defaultdict(Counter or Variant)"""
    if 'vars' in results[0]:  # variants are stored as defaultdict(Variant) indexed by Variant object
        master_var = defaultdict(kt.Variant)
        output = open((destination + results[0] + '.bed'), "a+")
        output.write("CHROM\tPOS\tREF\tALT\n")  # header same for all
        for variant in results[1:]:
            for k, v in variant.items():
                master_var[k] = v
                output.write(str(v))
        output.close()
        return master_var
    if 'transitions' in results[0]:  # kmer transitions stored as defaultdict(Counter) indexed by kmer
        master_count = defaultdict(Counter)
        for counts in results[1:]:
            for k, v in counts.items():
                for alt, count in v.items():
                    master_count[k][alt] += count
        merged_df = pd.DataFrame.from_dict(master_count, orient='index')
        fname = destination + results[0] + '.csv'
        merged_df.to_csv(fname)
        return master_count
    pass


class GenVCF:
    def __init__(self, ref_fasta_path, vcf_path, kmer_size, nprocs):
        self.vcf_path = vcf_path
        self.fasta_path = ref_fasta_path
        self.ref = Fasta(ref_fasta_path)
        self.vcf = VCF(vcf_path)
        self.kmer_size = kmer_size
        self.nprocs = nprocs
        self.keys = [c for c in self.vcf.seqnames if c in self.ref.keys()]
        self.directory = None
        if len(self.keys) == 0:
            self.keys = self.ref.keys()
            print('No common keys found. Using reference.')

    def set_destination(self, path):
        self.directory = path

    def full_string(self):
        ref_seq = ""
        for chrom in self.keys:
            ref_seq += str(self.ref[chrom])
        return ref_seq

    def get_kmer_frequency(self, klen):
        # returns a counter
        return kt.get_kmer_count(self.full_string(), klen)

    def vcf_scan(self):
        if self.nprocs == 0:
            self.nprocs = mp.cpu_count()
        regions = kt.get_split_vcf_regions(self.vcf_path, self.nprocs)
        args = [[region, self.vcf_path, self.fasta_path, self.kmer_size] for region in regions]
        pool = mp.Pool(self.nprocs)
        results = [funccall.get() for funccall in [pool.starmap_async(process_region, args)]]
        pool.close()
        all_vars, singletons, not_singletons = ['all_vars'], ['singleton_vars'], ['notsingleton_vars']
        all_transitions, singleton_transitions, notsingleton_transitions = ['all_transitions'], [
            'singleton_transitions'], ['notsingleton_transitions']
        for result in results[0]:
            for key, value in result.items():
                if key == 'all':
                    all_transitions.append(value[0])
                    all_vars.append(value[1])
                if key == 'singletons':
                    singleton_transitions.append(value[0])
                    singletons.append(value[1])
                if key == 'not_singletons':
                    notsingleton_transitions.append(value[0])
                    not_singletons.append(value[1])
        all_results = [all_vars, singletons, not_singletons, all_transitions, singleton_transitions,
                       notsingleton_transitions]
        destination = kt.prepare_directory(parent=self.directory)
        for res in all_results:
            merge_and_save(res, destination)
        return results


if __name__ == "__main__":
    start = time.time()
    f, v = kt.test_data()
    gv = GenVCF(f, v, 3, 6)
    r = gv.vcf_scan()
    print('Done in ' + str(time.time() - start) + ' seconds')
