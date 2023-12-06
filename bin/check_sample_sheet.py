#!/usr/bin/env python

from pathlib import Path
import logging

import pandas as pd
import typer
from rich.console import Console
from rich.logging import RichHandler


def main(input_path: Path,
         output_sample_sheet: Path):
    """Check and reformat sample sheet into CSV"""
    init_logging()
    logging.info(f'Input sample sheet: {input_path}')
    df = read_sample_sheet(input_path)
    fix_sample_names(df)
    fix_file_paths(df)
    df.to_csv(output_sample_sheet, index=False)
    logging.info(f'Wrote reformatted sample sheet CSV to "{output_sample_sheet}"')


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


def read_sample_sheet(input_path):
    ext = input_path.suffix.lower()
    logging.info(f'extension: {ext}')
    try:
        if ext in ['.tsv', '.txt', '.tab']:
            df = pd.read_table(input_path)
        elif ext == '.csv':
            df = pd.read_csv(input_path)
        elif ext in ['.xls', '.xlsx', '.ods']:
            df = pd.read_excel(input_path)
        else:
            raise ValueError(f'Unknown file format for sample sheet "{input_path}"')
    except Exception as ex:
        logging.exception(ex)
        raise ex
    logging.info(f'Sample sheet table shape: {df.shape}')
    assert df.shape[1] == 2, f'Two columns expected in sample sheet, but {df.shape[1]} found!'
    df.columns = ['sample', 'reads_path']
    return df


def fix_file_paths(df):
    reads_paths = []
    for fp in df.reads_path:
        if fp.startswith('http') or fp.startswith('ftp') or fp.startswith('s3'):
            reads_paths.append(fp)
        else:
            path = Path(fp)
            reads_paths.append(str(path.resolve().absolute()))
    df.reads_path = reads_paths


def fix_sample_names(df):
    samples: pd.Series = df['sample']
    contains_invalid_chars: pd.Series = samples.str.contains(r'[^\w\-]', regex=True)
    if contains_invalid_chars.any():
        linenos = contains_invalid_chars[contains_invalid_chars].index + 1
        logging.warning(f'{contains_invalid_chars.sum()} sample names contain invalid characters. Lines {list(linenos)} ')
        # remove leading and trailing non-word characters, replace remaining non-word characters with underscores
        samples = samples.str.replace(r'^[^\w\-]+', '', regex=True)
        samples = samples.str.replace(r'[^\w\-]+$', '', regex=True)
        samples = samples.str.replace(r'[^\w\-]+', '_', regex=True)
        for lineno, oldname, newname in zip(linenos, df['sample'][contains_invalid_chars], samples[contains_invalid_chars]):
            logging.warning(f'Changed sample name: [line {lineno}] "{oldname}" → "{newname}"')
        df['sample'] = samples


if __name__ == '__main__':
    typer.run(main)
