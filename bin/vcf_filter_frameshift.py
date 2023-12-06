#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging
from io import TextIOWrapper
from os import PathLike
from pathlib import Path
from typing import TextIO, Union, Tuple, Optional

import pandas as pd
import typer
from rich.console import Console
from rich.logging import RichHandler


VCF_COL_DTYPES: dict = dict(CHROM=str,
                            POS='uint32',
                            ID=str,
                            REF=str,
                            ALT=str,
                            QUAL=float,
                            FILTER=str,
                            INFO=str,
                            FORMAT=str)


def read_vcf(vcf_file: Union[PathLike, TextIOWrapper, TextIO]) -> Tuple[str, Optional[pd.DataFrame]]:
    """Read VCF file into a DataFrame"""
    header = ''
    with vcf_file if isinstance(vcf_file, (TextIOWrapper, TextIO)) else open(vcf_file) as fh:
        vcf_cols = []
        for line in fh:
            header += line
            if line.startswith('#CHROM'):
                vcf_cols = line[1:].strip().split('\t')
                break
        try:
            df = pd.read_table(fh,
                               comment='#',
                               header=None,
                               names=vcf_cols,
                               dtype=VCF_COL_DTYPES)
        except ValueError:
            df = None
    return header, df


def write_vcf(
    vcf_file: Union[PathLike, TextIOWrapper, TextIO], 
    header: str, 
    df: Optional[pd.DataFrame]
) -> None:
    with vcf_file if isinstance(vcf_file, (TextIOWrapper, TextIO)) else open(vcf_file, 'w') as fh:
        fh.write(header)
        if df is not None and not df.empty:
            df.to_csv(fh, sep='\t', header=False, index=False)


def main(input_vcf: Path = typer.Argument(None, help='VCF file to filter'), 
         output_vcf: Path = typer.Argument(None, help='Ouput VCF file to filter'),
         verbose: bool = typer.Option(True)):
    """Filter variants that lead to frameshift mutations from a VCF file

    Any indels resulting in a frameshift (indel length not divisible by 3 (AA codon length)) are
    removed from the output VCF file.
    """
    from rich.traceback import install

    install(show_locals=True, console=Console(stderr=True))

    logging.basicConfig(format='%(message)s',
                        datefmt='[%Y-%m-%d %X]',
                        level="NOTSET",
                        handlers=[RichHandler(rich_tracebacks=True,
                                              tracebacks_show_locals=True,
                                              console=Console(stderr=True))])
    log = logging.getLogger("rich")
    if input_vcf and not input_vcf.is_file():
        log.warning(f'input_vcf not a file or stdin stream!')
        sys.exit(1)

    header, df = read_vcf(input_vcf if input_vcf else sys.stdin)
    if df is None or df.empty:
        logging.warning(f'VCF has no variants! Either no reads mapped to the reference genome or the sample is identical to the reference.')
        df_filtered = None
    else:
        potential_frameshift_mask = df.apply(lambda x:  (len(x.ALT) - len(x.REF)) % 3 == 0, axis=1)
        df_filtered = df[potential_frameshift_mask]
        # log any potential frameshift variants that are going to be filtered out in output VCF
        if (~potential_frameshift_mask).any():
            log.info(f'{(~potential_frameshift_mask).sum()} frameshift variants filtered out')
            for _, row in df[~potential_frameshift_mask].iterrows():
                log.info(row)
    write_vcf(output_vcf if output_vcf else sys.stdout, header, df_filtered)


if __name__ == '__main__':
    typer.run(main)
