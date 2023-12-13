# CFIA-NCFAD/nf-virontus

[![GitHub Actions CI](https://github.com/CFIA-NCFAD/nf-virontus/actions/workflows/ci.yml/badge.svg)](https://github.com/CFIA-NCFAD/nf-virontus/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/CFIA-NCFAD/nf-virontus/actions/workflows/linting.yml/badge.svg)](https://github.com/CFIA-NCFAD/nf-virontus/actions/workflows/linting.yml)
[![Cite with Zenodo](https://zenodo.org/badge/229837060.svg)](https://zenodo.org/doi/10.5281/zenodo.10369022)

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

Oxford Nanopore viral sequence analysis pipeline.

## Introduction

This pipeline performs read mapping and variant calling with [Minimap2] and [Clair3]. A consensus sequence is generated from major variants and variants that would not cause potential frameshift mutations using [Bcftools] with masking of low coverage depth regions with `N` characters.

Optionally, amplicon primers can be trimmed with [iVar] if a BED file of primer coordinates is supplied.

If read mapping against the SARS-CoV-2 reference genome Wuhan-Hu-1 ([MN908947.3](https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3/)), [Pangolin] lineage assignment will be done, as well as, [Nextclade] analysis.

## Pipeline Overview

>*NOTE:* This pipeline is still a work-in-progress. The following diagram shows the planned features and steps of the workflow:

```mermaid
flowchart LR
    classDef input fill:#fcba64,color:black
    classDef output fill:#b1fc9c,color:black
    subgraph legend["<b>Legend"]
        style legend fill:white,fill-opacity:0.5
        input([Input]):::input
        output([Output]):::output
        process[Process]
    end
    
    subgraph rqc["<b>fa:fa-dna Reads QA/QC</b>"]
        A([fa:fa-file Filtered Reads FASTQ]):::output
        FR["fa:fa-check Read QC & Filtering <br><small> fastp, nanoqc"]
        HR["fa:fa-cancel Dehosting <br><small> Kraken2 (optional)"]
        RR([fa:fa-file Raw Reads FASTQ]):::input --> FR 
        FR --> HR
        HR --> A
    end
    rqc --> rs
    subgraph rs[<b>fa:fa-crosshairs Reference Selection]
        RS["fa:fa-filter Ref Seq Selection <br><small> de novo assembly & BLAST, Mash, Kraken2"]
        frs_rs([fa:fa-file Reads FASTQ]):::input
        frs_rs --> RS
        R([fa:fa-file Ref Seqs FASTA]):::input --> RS
        RS --> T([fa:fa-file Top Ref Seq FASTA]):::output
    end

    rs --> rma
    rqc --> rma

    subgraph rma[<b>fa:fa-industry Reference Mapped Assembly]
        direction TB
        trs([fa:fa-file Top Ref FASTA]):::input
        fr([fa:fa-file FASTQ]):::input
        RM["fa:fa-bars-staggered Read Mapping <br><small> Minimap2"]
        PT["fa:fa-scissors Primer Trimming <br><small> iVar (optional)"]
        VC["fa:fa-code-compare Variant Calling <br><small> Clair3, Medaka"]
        VE["fa:fa-flask Variant Effect <br><small> SnpEff, SnpSift (optional)"]
        BAM([fa:fa-file BAM]):::output
        D["fa:fa-chart-area Coverage Stats <br><small> Mosdepth, Samtools"]
        CS["fa:fa-code-merge Make Consensus Sequence <br><small> Bcftools"]
        vcf([fa:fa-file VCF]):::output
        muts([fa:fa-table Amino Acid Mutations]):::output
        covbed([fa:fa-file Coverage BED]):::output
        fr --> RM
        trs --> RM
        trs --> VC
        RM --> BAM
        BAM --> PT
        PT --> BAM
        BAM --> D
        BAM --> VC
        vcf --> VE
        VE --> vcf
        VE --> muts
        VC --> vcf
        vcf --> CS
        trs --> CS
        D --> covbed
        covbed --> CS
        CS --> csf([fa:fa-file Consensus Sequence FASTA]):::output
    end
    rqc --> reporting
    rs --> reporting
    rma --> reporting
    subgraph reporting[<b>fa:fa-clipboard Reporting & Visualization]
        MQC[fa:fa-stethoscope MultiQC]
        MQCR([fa:fa-file MultiQC HTML Report]):::output
        CP[fa:fa-chart-area Seq Coverage Plots]
        png([fa:fa-image PNG]):::output
        pdf([fa:fa-file PDF]):::output
        MQC --> MQCR
        CP --> png & pdf
    end
```

## Installation

You will need to install [Nextflow] in order to run the Virontus pipeline.

> **NB:** [Singularity] or [Docker] is recommended for portable and reproducible execution of the pipeline with the `-profile singularity` or `-profile docker` command-line argument.

### 1) Install [Nextflow]

If you have [Conda] installed, you can install [Nextflow] with the following command:

```bash
conda install -c bioconda -c conda-forge nextflow
```

### 2) Install [Docker][] and/or [Singularity][]

Installing [Docker][] and/or [Singularity] is optional but recommended for portability and reproducibility of results.

### 3) Install Virontus

Nextflow will automatically download the latest version of Virontus. You can show the Virontus help message with usage information with:

```bash
nextflow run CFIA-NCFAD/nf-virontus --help
```

## Usage

Basic usage for mapping to SARS-CoV-2 reference genome [MN908947.3](https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3/) and [ARTIC V3 protocol](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V3) primers:

```bash
nextflow run CFIA-NCFAD/nf-virontus \
  --input samplesheet.csv \
  --genome MN908947.3 \
  --bed artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.bed
```

Can be simplified with:

```bash
nextflow run CFIA-NCFAD/nf-virontus \
  --input samplesheet.csv \
  --scov2 \
  --artic_v3
  # or `--freed` for Freed et al (2020) 1200bp amplicon method
```

Show usage information with

```bash
nextflow run CFIA-NCFAD/nf-virontus --help
```

> **NB:** See the [usage docs](docs/usage.md) for more info.

## Output

See the [output docs](docs/output.md) for more info.

## Credits

CFIA-NCFAD/nf-virontus was originally written by Peter Kruczkiewicz.

Bootstrapped with [nf-core/tools](https://github.com/nf-core/tools) `nf-core create`.

Thank you to the [nf-core/tools](https://github.com/nf-core/tools) team for a great tool for bootstrapping creation of a production ready Nextflow workflows.

[Bcftools]: https://samtools.github.io/bcftools/bcftools.html
[Conda]: https://conda.io/
[Docker]: https://www.docker.com/
[iVar]: https://github.com/andersen-lab/ivar
[Minimap2]: https://github.com/lh3/minimap2
[Nextflow]: https://www.nextflow.io
[Nextclade]: https://clades.nextstrain.org/
[Pangolin]: https://github.com/cov-lineages/pangolin/
[Singularity]: https://github.com/apptainer/apptainer
