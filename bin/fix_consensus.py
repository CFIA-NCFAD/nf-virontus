#!/usr/bin/env python
# coding: utf-8
import click
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


@click.command()
@click.option('-s', '--consensus-fasta', required=True,
              help='Consensus sequence not accounting for low coverage positions.')
@click.option('-d', '--depths-table', required=True,
              help='Samtools depth output table containing fields: ["barcode", "genome", "position", "depth"].')
@click.option('-o', '--output-fasta', required=True,
              help='Top mapping reference sequence FASTA output path.')
@click.option('--cov-threshold', default=1, help='Coverage depth threshold.')
@click.option('--low-cov-char', default='N', help='Low coverage depth position character replacement (default: "N")')
@click.option('--sample-name', default='', help='Optional sample name for FASTA header')
def fix(consensus_fasta, depths_table, output_fasta, cov_threshold, low_cov_char, sample_name):
    """Replace low coverage positions in a consensus sequence with N or some other character.

    The consensus sequence from `bcftools consensus` may have positions that 
    are incorrectly called as the reference when they should be called more 
    ambiguously as "N" or "-".
    """
    seq_rec = SeqIO.read(consensus_fasta, format='fasta')
    df = pd.read_csv(depths_table, 
                     sep='\t', 
                     header=None, 
                     names=['barcode', 
                            'genome', 
                            'position', 
                            'depth'])
    mutable_seq = seq_rec.seq.tomutable()
    for position in df[df.depth < cov_threshold].position:
        mutable_seq[position - 1] = low_cov_char

    consensus_seq_rec = SeqRecord(id=sample_name or seq_rec.id,
                                  seq=mutable_seq.toseq(),
                                  description=f'consensus_from_ref="{seq_rec.id}"')
    SeqIO.write([consensus_seq_rec], output_fasta, 'fasta')


if __name__ == '__main__':
    fix()
