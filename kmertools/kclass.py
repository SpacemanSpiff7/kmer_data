from collections import defaultdict, Counter


class Variant:
    def __init__(self, *args, **kwargs):
        if len(kwargs) == 0:
            self.REF = args[0]
            self.ALT = args[1]
            self.POS = args[2]
            self.CHROM = args[3]
        else:
            variant = kwargs.get('variant')
            self.REF = variant.REF
            self.ALT = variant.ALT
            self.POS = variant.POS
            self.CHROM = variant.CHROM

    def __str__(self):
        # return "POS: " + str(self.POS) + "\tREF: " + str(self.REF) + "\tALT: " + str(self.ALT) + '\n'
        return str(self.CHROM) + "\t" + str(self.POS) + "\t" + str(self.REF) + "\t" + str(self.ALT) + '\n'

    def __repr__(self):
        return "CHROM: " + str(self.CHROM) + "\t" + "POS: " + str(self.POS) + "\tREF: " + str(
            self.REF) + "\tALT: " + str(self.ALT) + '\n'

    def __eq__(self, other):
        if not isinstance(other, Variant):
            return False
        return other.CHROM == self.CHROM and other.POS == self.POS and other.REF == self.REF and other.ALT == self.ALT

    def __ne__(self, other):
        if not isinstance(other, Variant):
            return True
        return not (
                other.CHROM == self.CHROM and other.POS == self.POS and other.REF == self.REF and other.ALT == self.ALT)

    def __hash__(self):
        return hash(str(self))

    def __lt__(self, other):
        if other.CHROM != self.CHROM:
            return self.CHROM < other.CHROM
        else:
            return self.POS < other.POS

    def __gt__(self, other):
        if other.CHROM != self.CHROM:
            return self.CHROM > other.CHROM
        else:
            return self.POS > other.POS

    def csv_str(self):
        return str(self.CHROM) + "\t" + str(self.POS) + "\t" + str(self.REF) + "\t" + str(self.ALT) + '\n'


class Kmer:
    def __init__(self, sequence):
        from kmertools.kmethods import get_complementary_sequence
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


class VCFRegion:
    def __init__(self, chrom, start, stop):
        self.region = [chrom, start, stop]

    # def add(self, chrom, start, stop):
    #     if chrom in self.regions.keys():
    #         self.regions[chrom][1] = stop

    def size(self):
        return self.region[2] - self.region[1]

    def __str__(self):
        return str(self.region[0]) + ':' + str(self.region[1]) + '-' + str(self.region[2])

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(str(self))


class Storage:
    def __init__(self, *args, **kwargs):
        self.items = []
        self.name = args[1]
        for i in range(2, len(args)):
            self.items.append(i)

    def __hash__(self):
        return hash(self.name + str(len(self.items)))

    def __lt__(self, other):
        return hash(self) - hash(other)

    def __eq__(self, other):
        return hash(self) == hash(other)

    def get_positions(self):
        for i in self.items:
            if isinstance(i, defaultdict(Variant)):
                return i
        return False

    def get_transitions(self):
        for i in self.items:
            if isinstance(i, defaultdict(Counter)):
                return i
        return False
