#!/usr/bin/env python3

import os
import shutil
import sys

def fix_bad_header(filep):
    from_file = open(filep, 'r')
    line = from_file.readline()
    to_file = open(filep, 'w')
    nline = line.split('REF_vcf')
    to_file.write(line[:-len(nline[1])] + '\n')
    to_file.write(nline[1])
    shutil.copyfileobj(from_file, to_file)
    from_file.close()
    to_file.close()

if __name__=="__main__":
    fp = sys.argv[1]
    fix_bad_header(fp)
