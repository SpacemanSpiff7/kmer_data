from collections import defaultdict
import multiprocessing as mp
from cyvcf2 import VCF
from kmertools import is_quality_singleton, generate_csv_from_variants
from kmertools.kclass import Variant

vcf_path = '../samp_build38.vcf.bgz'


def process_region(chr_num):
    vcf = VCF(vcf_path)
    variant_positions = defaultdict(Variant)
    for variant in vcf(chr_num):
        if is_quality_singleton(variant):
            variant_positions[variant.POS] = Variant(variant.REF, "".join(variant.ALT), variant.POS, variant.CHROM)
    print("Done processing chromosome " + str(chr_num))
    return variant_positions


# def run_process(nthreads):
#     chromosomes = list(REF_GENOME.keys())
#     pool = mp.Pool(nthreads)
#     results = [funccall.get() for funccall in [pool.map_async(process_region, chromosomes)]]
#     pool.close()
#     return results


# def parallel_VCF_read(vcf_fpath, nthreads=0, outfile='variants_samp.csv'):
#     # TODO: split up VCF into equal parts
#     if nthreads == 0:
#         nthreads = mp.cpu_count()
#     global vcf_path
#     vcf_path = vcf_fpath
#     results = run_process(nthreads)
#     output = open(outfile, "a+")
#     output.write("CHROM\tPOS\tREF\tALT\n")  # header same for all
#     for dic in results[0]:
#         for k, v in dic.items():
#             output.write(str(v))
#     output.close()
#     print("Done reading VCF file.")
#
# # if __name__ == "__main__":
# #     parallel_variant()
