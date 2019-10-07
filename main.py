from kmertools import *
import sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--numthreads', '-N', action='store', dest='N_THREADS', help='Number of threads available',
                        default=4)
    parser.add_argument('--filepath', '-f', action='store', dest='vcf_path',
                        help='Enter path to VCF file or CSV file containing variants', default=None)
    parser.add_argument('--kmer-len', '-k', action='store', dest='kmer_length', help='Enter length of Kmer', default=3)
    parser.add_argument('--output', '-o', action='store', dest='output_file', help='Enter desired output file name',
                        default='variants_samp.csv')
    args = parser.parse_args()
    try:
        nthreads = int(args.N_THREADS)
        vcf_path = args.vcf_path
        kmer_len = int(args.kmer_length)
        o_file = args.output_file
    except ValueError:
        print("Invalid parameters. Exiting...")
        exit(0)
    parallel_variant(vcf_path, nthreads=nthreads, outfile=o_file)
