#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from pathlib import Path

import typer
import pandas as pd
from rich.console import Console
from rich.logging import RichHandler


def init_logging():
    from rich.traceback import install
    install(show_locals=True,
            width=120,
            word_wrap=True,
            console=Console(stderr=True))
    logging.basicConfig(format='%(message)s',
                        datefmt='[%Y-%m-%d %X]',
                        level=logging.DEBUG,
                        handlers=[RichHandler(rich_tracebacks=True,
                                              tracebacks_show_locals=True,
                                              console=Console(stderr=True))])


def main(
    input_bed: Path,
    output_bed: Path,
    left_suffix: str = typer.Option('_LEFT', help='Left primer suffix'),
    right_suffix: str = typer.Option('_RIGHT', help='Right primer suffix')
) -> None:
    init_logging()
    logging.info(f'Reading "{input_bed}" into dataframe')
    df_primer = pd.read_table(input_bed, 
        header=None, 
        names=['genome', 'start_idx', 'end_idx', 'primer', 'pool', 'strand'])
    logging.info(f'Removing left/right primer suffixes "{left_suffix}" and "{right_suffix}" from primer IDs.')
    df_primer['amplicon'] = df_primer.primer.str.replace(f'({left_suffix}|{right_suffix}).*$', '', regex=True)
    logging.info(f'Grouping primers by amplicon names and building amplicon BED file.')
    df_amplicon = df_primer.groupby('amplicon').agg(
        genome=('genome', 'first'),
        start=('start_idx', 'min'),
        end=('end_idx', 'max'),
        amplicon=('amplicon', 'first'),
        pool=('pool', 'first')
    )
    df_amplicon['strand'] = '+'
    df_amplicon.sort_values('start', inplace=True)
    logging.info(f'Writing amplicon BED file to "{output_bed}"')
    df_amplicon.to_csv(output_bed, index=False, header=None, sep='\t')
    logging.info('Done!')


if __name__ == '__main__':
    typer.run(main)
