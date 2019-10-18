REF_GENOME = dict({'1': 249250621,
                   '10': 135534747,
                   '11': 135006516,
                   '12': 133851895,
                   '13': 115169878,
                   '14': 107349540,
                   '15': 102531392,
                   '16': 90354753,
                   '17': 81195210,
                   '18': 78077248,
                   '19': 59128983,
                   '2': 243199373,
                   '20': 63025520,
                   '21': 48129895,
                   '22': 51304566,
                   '3': 198022430,
                   '4': 191154276,
                   '5': 180915260,
                   '6': 171115067,
                   '7': 159138663,
                   '8': 146364022,
                   '9': 141213431,
                   'X': 155270560,
                   'Y': 59373566})

REF_FASTA_PATH = "/Users/simonelongo/too_big_for_icloud/ref_genome/REFERENCE_GENOME_GRch37.fa"

VCF_PATH = "/Users/simonelongo/Documents/QuinlanLabFiles/kmer_data/samp.vcf.bgz"

KMER_SIZE = 3


def set_kmer_size(length):
    global KMER_SIZE
    KMER_SIZE = length
    return


def set_vcf_path(path):
    global VCF_PATH
    VCF_PATH = path
    return
