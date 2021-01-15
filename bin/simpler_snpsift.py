#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging
from io import TextIOWrapper
from os import PathLike
from pathlib import Path
from typing import TextIO, Union, Tuple

import pandas as pd
import typer
from rich.console import Console
from rich.logging import RichHandler


def main(snpsift_input: Path = typer.Argument(None, help='SnpSift tab-delimited table of variant effects'), 
         output: Path = typer.Argument(None, help='Ouput simple SnpSift table')):
    df = pd.read_table(snpsift_input)
    series = []
    for c in df.columns:
        if c == 'AC':
            df_ac = df[c].str.split(',', n=1, expand=True)
            REF_AC = df_ac[0].astype(int)
            REF_AC.name = 'REF_AC'
            ALT_AC = df_ac[1].astype(int)
            ALT_AC.name = 'ALT_AC'
            AF = ALT_AC / (REF_AC + ALT_AC)
            AF.name = 'AF'
            series += [REF_AC, ALT_AC, AF]
            continue
        idx = c.find('[*].')
        if idx > 0:
            new_series = df[c].str.split(',', n=1, expand=True)[0]
            new_series.name = c[idx+4:].lower()
            series.append(new_series)
        else:
            series.append(df[c])
    df_out = pd.concat(series, axis=1)
    df_out.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    typer.run(main)