from collections import defaultdict
import psutil
import multiprocessing as mp
from cyvcf2 import VCF
from kmertools import REF_GENOME, is_quality_variant, generate_csv_from_variants
from kmertools.kclass import Variant


def process_region(chr_num):
    vcf_path = '../samp.vcf.bgz'
    vcf = VCF(vcf_path)
    variant_positions = defaultdict(Variant)
    for variant in vcf(chr_num):
        if is_quality_variant(variant):
            variant_positions[variant.POS] = Variant(variant.REF, "".join(variant.ALT), variant.POS, variant.CHROM)
    return variant_positions


def parallel_variant(nthreads=0, outfile='variants_samp.csv'):
    if nthreads == 0:
        nthreads = psutil.cpu_count(logical=False)
    chromosomes = list(REF_GENOME.keys())
    pool = mp.Pool(nthreads)
    results = [funccall.get() for funccall in [pool.map_async(process_region, chromosomes)]]
    pool.close()
    output = open(outfile, "a+")
    output.write("CHROM\tPOS\tREF\tALT\n")  # header same for all
    for dic in results[0]:
        for k, v in dic.items():
            output.write(str(v))
    output.close()


if __name__ == "__main__":
    parallel_variant()