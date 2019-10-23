from pyfaidx import Fasta
from cyvcf2 import VCF
import kmertools as kt


class VCFScanner:
    def __init__(self, ref_fasta, vcf):
        self.ref = Fasta(ref_fasta)
        self.vcf = VCF(vcf)
        self.names = kt.get_fasta_primary_chroms(ref_fasta)

    def variant_kmer(self, kmer_size, singleton=True, non_sing=True, keys=None):
        if keys is None:
            keys = self.names

