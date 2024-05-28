# CFIA-NCFAD/nf-flu changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[2.0.2dev](https://github.com/CFIA-NCFAD/nf-virontus/releases/tag/2.0.2dev)]

* hotfix for properly parsing primer BED file with extra `sequence` column

## [[2.0.1dev](https://github.com/CFIA-NCFAD/nf-virontus/releases/tag/2.0.1dev)]

This release fixes an issue where random ambiguous bases are introduced into the consensus sequence by Bcftools (#9)

## [[2.0.0dev](https://github.com/CFIA-NCFAD/nf-virontus/releases/tag/2.0.0dev)] - 2023-12-06

This release is a major refactor using Nextflow DSL-2 modules and similar process structure and definitions to nf-core pipelines. A single reference genome is assumed for all samples unlike previous versions of nf-virontus that would perform parallel analyses on all provided reference sequences. A smart reference genome selection module is in the works.

For SARS-Cov-2 analysis, the `--scov2` flag can be used to specify the MN908947.3 Wuhan-Hu-1 sequence as the reference genome. This is equivalent to specifying `--genome MN908947.3`. Pangolin and Nextclade analysis will automatically be performed if the reference genome is MN908947.3.

**NOTE:** This is a in-development release that does not have all features expected for version 2.0.0, but should provide output that is compatible with the [xlavir](https://github.com/CFIA-NCFAD/xlavir) Excel Nextflow bioinformatics pipeline report generation tool.

### Changes

* refactor into DSL-2 modules
* Clair3 variant calling instead of Medaka/Longshot for faster more accurate variant calling esp. at lower coverage depths
* adapted nf-core/modules processes (Bcftools, Nextclade, Nextalign, PycoQC, NanoPlot, MultiQC)
* added SnpEff variant effect prediction
* added SnpSift to parse SnpEff annotations into tab-delimited table like nf-core/viralrecon
* allele fraction threshold set at 0.75 by default for major variant to be incorporated into consensus sequence
* added `bin/simpler_snpsift.py` to simplify SnpEff/SnpSift table
* updated tools and dependencies to latest versions

## 1.1.0

* Medaka and Longshot for variant calling.

> *NOTE:* This version is no longer supported and is not recommended for analyses.
