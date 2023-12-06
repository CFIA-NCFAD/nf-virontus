#!/usr/bin/env python

import gzip
from pathlib import Path
import sys
import re


def to_vcf_info(form: str, sample: str) -> str:
    """Convert Clair3 FORMAT and SAMPLE columns to VCF INFO column

    >>> to_vcf_info('GT:DP:AD:AF', '0/1:3:0,3:1.0')
    'GT=0/1;DP=3;AD=0,3;AF=1.0'
    """
    return ';'.join(f'{f}={s}' for f, s in zip(form.split(':'), sample.split(':')))


def main():
    """Add DP,AF info to INFO column in Clair3 VCF

    The DP and AF information is found in the FORMAT and SAMPLE columns (last 2 columns). To make the VCF compatible
    with SnpEff/SnpSift, it's necessary to add this info to the INFO column.

    Output is written to stdout.
    """
    vcf_gz_path = Path(sys.argv[1])
    with gzip.open(vcf_gz_path) as f:
        for l in f:
            line = l.decode().strip()
            # skip lines containing Clair3 RefCall using regex search
            if re.search(r'(?<=\t)RefCall(?=\t)', line):
                continue
            # skip lines starting with '#' (header lines)
            if line.startswith('#'):
                print(line)
            else:
                sp = line.split('\t')
                *_, info, format, sample = sp
                # add DP,AF info to INFO column
                sp[7] = f'{info};{to_vcf_info(format, sample)}'
                print('\t'.join(sp))


if __name__ == '__main__':
    main()
