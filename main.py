import kmertools as km
from kmertools.workflows import vcf_process
import sys
import argparse
import time


# /scratch/general/lustre/u0319040/ref_genome/gnomad.genomes.r2.1.1.sites.vcf.bgz
def main():
    # TODO: Implement ref directory check
    if fasta_path is not None:
        km.REF_FASTA_PATH = fasta_path
    if vcf_path is not None:
        km.VCF_PATH = vcf_path
        if kmer_len > 0 and kmer_len != km.constants.KMER_SIZE:
            km.KMER_SIZE = kmer_len
        # km.parallel_VCF_read(vcf_path, nthreads=nthreads, outfile=o_file)
        km.vcf_parallel(vcf_path, nthreads, outfile=o_file)
        var_path = o_file
        # if kmer_len > 0:
        #     if kmer_len != km.constants.KMER_SIZE:
        #         km.KMER_SIZE = kmer_len
        #     if vcf_path is None and csv_fpath is not None:
        #         var_path = csv_fpath  # CSV is then read appropriately in 'process variants'
        #     else:
        #         print('No input file to work from, please reevaluate arguments.')
        #         exit(1)
        #     if kmer_len != km.constants.KMER_SIZE:
        #         km.KMER_SIZE = kmer_len
        #     # var_counts = km.process_variants_old(var_path, kmer_len, fasta_path)  # Pandas data frame
        #     # km.find_ref_kmer_freq(kmer_len, fasta_path)
    return True


if __name__ == "__main__":
    start = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('--singleton-flow', '-fs', action='store', dest='workflow_singleton', default=None)
    parser.add_argument('--numthreads', '-N', action='store', dest='N_THREADS', help='Number of threads available',
                        default=4)
    parser.add_argument('--ref-directory', '-rd', action='store', dest='ref_dir_path', default=None,
                        help='Enter path to directory containing indexed VCF and Fasta files.')
    parser.add_argument('--vcf-path', '-v', action='store', dest='vcf_path',
                        help='Enter path to VCF file or CSV file containing variants', default=None)
    parser.add_argument('--variant-csv', '-t', action='store', dest='csv_var_fpath',
                        help='Enter filepath for CSV file containing variants if available', default=None)
    parser.add_argument('--reference-genome', '-rg', action='store', dest='ref_fasta',
                        help='Enter path to reference genome fasta file',
                        default='/Users/simonelongo/too_big_for_icloud/REFERENCE_GENOME_GRch37.fa')
    parser.add_argument('--kmer-len', '-k', action='store', dest='kmer_length', help='Enter length of Kmer', default=3)
    parser.add_argument('--output', '-o', action='store', dest='output_file', help='Enter desired output file name',
                        default='variants_samp.csv')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    try:
        is_wkflow = args.workflow_singleton is not None
        nthreads = int(args.N_THREADS)
        vcf_path = args.vcf_path
        kmer_len = int(args.kmer_length)
        o_file = args.output_file
        csv_fpath = args.csv_var_fpath
        fasta_path = args.ref_fasta
        ref_dir = args.ref_dir_path
    except ValueError:
        print("Invalid parameters. Exiting...")
        exit(1)
    if is_wkflow:
        vcf_process.vcf_singleton_analysis(vcf_path, fasta_path, kmer_len, nprocs=nthreads, var_out_file=o_file)
        print("Done in " + str(time.time() - start))
        exit(0)
    main()
