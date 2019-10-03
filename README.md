# Kmer Genome Search Tools
 Written by Simone Longo October 3, 2019 
at the University of Utah
----

This package contains tools for manipulating VCF files building on the work done by **cyvcf2**

It also contains convenient classes to use when considering DNA sequences.

```
from cyvcf2 import VCF
from kclass import Variant
```
`class Variant` is initialized from a variant object returned bya cyvcf2.VCF iterator or from entering the values for ALT, REF, POS, and CHROM in that order
For example: 

`variant = Variant('A', 'T', 123456, 19)`

OR:

```
for variant in VCF(vcf_filepath):
    var_object = Variant(variant)
```

There is also a class for kmer sequences, `Kmer`

`from kclass import Kmer`

This class treats sequences and their reverse complements as identical and stores both in the same object

```
k1 = Kmer('ACTG')
k2 = Kmer('CAGT')
k1 == k2
```
```
Out[3]: True
```


