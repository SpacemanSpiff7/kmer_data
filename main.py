import kmertools as km
import sys
import argparse


def main():
    if vcf_path is not None:
        km.parallel_variant(vcf_path, nthreads=nthreads, outfile=o_file)
        var_path = o_file
    if kmer_len > 0:
        if vcf_path is None and csv_fpath is not None:
            var_path = csv_fpath
        else:
            print('No input file to work from, please reevaluate arguments.')
            exit(1)
        var_counts = km.process_variants(var_path, kmer_len)  # Pandas dataframe
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--numthreads', '-N', action='store', dest='N_THREADS', help='Number of threads available',
                        default=4)
    parser.add_argument('--vcf-filepath', '-vf', action='store', dest='vcf_path',
                        help='Enter path to VCF file or CSV file containing variants', default=None)
    parser.add_argument('--variant-csv', '-cf', action='store', dest='csv_var_fpath',
                        help='Enter filepath for CSV file containing variants if available', default=None)
    parser.add_argument('--kmer-len', '-k', action='store', dest='kmer_length', help='Enter length of Kmer', default=0)
    parser.add_argument('--output', '-o', action='store', dest='output_file', help='Enter desired output file name',
                        default='variants_samp.csv')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    try:
        nthreads = int(args.N_THREADS)
        vcf_path = args.vcf_path
        kmer_len = int(args.kmer_length)
        o_file = args.output_file
        csv_fpath = args.csv_var_fpath
    except ValueError:
        print("Invalid parameters. Exiting...")
        exit(1)
    main()
