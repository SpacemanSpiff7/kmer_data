from src.kmethods import get_complementary_sequence


class Variant:
    def __init__(self, ref, alt, pos, chrom=22):
        self.REF = ref
        self.ALT = alt
        self.POS = pos
        self.chrom = chrom

    def __str__(self):
        # return "POS: " + str(self.POS) + "\tREF: " + str(self.REF) + "\tALT: " + str(self.ALT) + '\n'
        return str(self.chrom) + "\t" + str(self.POS) + "\t" + str(self.REF) + "\t" + str(self.ALT) + '\n'

    def __repr__(self):
        return "CHROM: " + str(self.chrom) + "\t" + "POS: " + str(self.POS) + "\tREF: " + str(
            self.REF) + "\tALT: " + str(self.ALT) + '\n'

    def csv_str(self):
        return str(self.chrom) + "\t" + str(self.POS) + "\t" + str(self.REF) + "\t" + str(self.ALT) + '\n'


class Kmer:
    def __init__(self, sequence):
        self.sequence = sequence
        self.complement = get_complementary_sequence(sequence)
        if self.complement < self.sequence:
            temp = self.sequence
            self.sequence = self.complement
            self.complement = temp

    def __eq__(self, other):
        return other.sequence == self.sequence or other.sequence == self.complement

    def __ne__(self, other):
        return other.sequence != self.sequence and other.sequence != self.complement

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return self.sequence

    def __hash__(self):
        return hash(self.sequence)

    def __lt__(self, other):
        return other.sequence < self.sequence

    def __gt__(self, other):
        return other.sequence > self.sequence
