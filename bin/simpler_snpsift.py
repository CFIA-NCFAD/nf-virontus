#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import re
from pathlib import Path

import pandas as pd
import typer
from rich.console import Console
from rich.logging import RichHandler

aa_codes = dict(
    ALA='A',
    ARG='R',
    ASN='N',
    ASP='D',
    CYS='C',
    GLU='E',
    GLN='Q',
    GLY='G',
    HIS='H',
    ILE='I',
    LEU='L',
    LYS='K',
    MET='M',
    PHE='F',
    PRO='P',
    SER='S',
    THR='T',
    TRP='W',
    TYR='Y',
    VAL='V',
)


def parse_aa(
    gene: str,
    ref: str,
    alt: str,
    nt_pos: int,
    aa_pos: int,
    snpeff_aa: str,
    effect: str
) -> str:
    if snpeff_aa == '.':
        return f'{ref}{nt_pos}{alt}'
    m = re.match(r'p\.([a-zA-Z]+)(\d+)([a-zA-Z]+)', snpeff_aa)
    if m is None or snpeff_aa.startswith('p.'):
        aa_str = snpeff_aa[2:].upper()
        for aa3, aa1 in aa_codes.items():
            aa_str = aa_str.replace(aa3, aa1)
        if aa_str.endswith('DEL'):
            aa_str = aa_str.replace('DEL', 'del')
        return f'{ref}{nt_pos}{alt}({gene}:{aa_str})'
    ref_aa, aa_pos_str, alt_aa = m.groups()
    ref_aa = get_aa(ref_aa)

    if effect == 'stop_lost':
        alt_aa = get_aa(alt_aa.replace('ext', ''))
        return f'{ref}{nt_pos}{alt}({gene}:{ref_aa}{aa_pos_str}{alt_aa}[stop_lost])'
    if effect == 'frameshift_variant':
        return f'{ref}{nt_pos}{alt}({gene}:{ref_aa}{aa_pos_str}[FRAMESHIFT])'
    if effect == 'conservative_inframe_deletion':
        return f'{ref}{nt_pos}{alt}({gene}:{ref_aa}{aa_pos_str}{alt_aa})'
    if effect == 'disruptive_inframe_deletion':
        return f'{ref}{nt_pos}{alt}({gene}:{ref_aa}{aa_pos_str}{alt_aa}[disruptive_inframe_deletion])'

    alt_aa = get_aa(alt_aa)
    return f'{ref}{nt_pos}{alt}({gene}:{ref_aa}{aa_pos_str}{alt_aa})'


def get_aa(s: str) -> str:
    out = ''
    for i in range(0, len(s), 3):
        aa = s[i: i + 3]
        try:
            aa_code = aa_codes[aa.upper()]
        except KeyError:
            aa_code = aa
        out += aa_code
    return out


def main(
        snpsift_input: Path = typer.Argument(None, help='SnpSift tab-delimited table of variant effects'),
        output: Path = typer.Argument(None, help='Ouput simple SnpSift table')
    ):
    """Simplify the SnpSift output from annotated Medaka and SnpEff VCF

    SR INFO tag contains "Depth of spanning reads by strand which best align to 
    each allele (ref fwd, ref rev, alt1 fwd, alt1 rev, etc.)".
    This info is used to calculate the REF allele count (AC) and ALT AC.
    AF is calculated as ALT_AC / (REF_AC + ALT_AC).
    Ambiguous bases/alleles are ignored in the AF calculation.
    """
    from rich.traceback import install
    install(show_locals=True,
            width=120,
            word_wrap=True,
            console=Console(stderr=True))
    logging.basicConfig(
        format='%(message)s',
        datefmt='[%Y-%m-%d %X]',
        level=logging.DEBUG,
        handlers=[
            RichHandler(
                rich_tracebacks=True,
                tracebacks_show_locals=True,
                console=Console(stderr=True),
            )
        ]
    )

    df = pd.read_table(snpsift_input)
    if not df.empty:
        series = []
        for c in df.columns:
            if c == 'SR':
                df_ac = df[c].str.split(',', expand=True)
                REF_AC = df_ac[0].astype(int) + df_ac[1].astype(int)
                REF_AC.name = 'REF_AC'
                ALT_AC = df_ac[2].astype(int) + df_ac[3].astype(int)
                ALT_AC.name = 'ALT_AC'
                AF = ALT_AC / (REF_AC + ALT_AC)
                AF.name = 'AF'
                series += [REF_AC, ALT_AC, AF]
                continue
            idx = c.find('[*].')
            if idx > 0:
                new_series = df[c].astype(str).str.split(',', n=1, expand=True)[0]
                new_series.name = c[idx + 4:].lower()
                series.append(new_series)
            else:
                series.append(df[c])
        df_out = pd.concat(series, axis=1)
        mutation_desc = []
        for row in df_out.itertuples():
            mutation_desc.append(
                parse_aa(
                    gene=row.gene,
                    ref=row.REF,
                    alt=row.ALT,
                    nt_pos=row.POS,
                    aa_pos=row.aa_pos,
                    snpeff_aa=row.aa,
                    effect=row.effect
                )
            )
        df_out['mutation'] = mutation_desc
        df_out.to_csv(output, sep='\t', index=False)
    else:
        expected_columns = [
            'CHROM',
            'POS',
            'REF',
            'ALT',
            'DP',
            'REF_AC',
            'ALT_AC',
            'AF',
            'gene',
            'geneid',
            'impact',
            'effect',
            'codon',
            'aa',
            'cdna_pos',
            'cdna_len',
            'cds_pos',
            'cds_len',
            'aa_pos',
            'aa_len',
            'mutation',
        ]
        with open(output, 'w') as fout:
            expected_header = "\t".join(expected_columns)
            fout.write(f'{expected_header}\n')


if __name__ == '__main__':
    typer.run(main)
