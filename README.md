# peterk87/nf-virontus
**Virontus viral Oxford Nanopore sequence analysis pipeline**

[![Build Status](https://travis-ci.org/peterk87/nf-virontus.svg?branch=master)](https://travis-ci.org/peterk87/nf-virontus)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/peterk87/nf-virontus.svg)](https://hub.docker.com/r/peterk87/nf-virontus)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/4297)


### Introduction

The Virontus pipeline is for the analysis of viral shotgun and amplicon  Oxford Nanopore sequence data. Given basecalled (and demultiplexed) Nanopore reads, Virontus produces one or more consensus sequences from read mapping with [Minimap2] and variant calling with [Medaka] and [Longshot] results with respect to one or more reference sequences. For amplicon sequencing, the user should provide a BED file containing primer coordinates with respect to a reference sequence so that the primer sequences can be trimmed using [iVar]. 


Optionally, Virontus will perform taxonomic classification with [Kraken2] and [Centrifuge] if index paths are provided. Reads can be filtered by taxonomic classification. By default viral and unclassified reads are filtered.

De novo assembly with [Unicycler] can be optionally performed if desired (specify `--do_unicycler_assembly` when running Virontus). If taxonomic classification is performed then taxonomically filtered reads will be assembled, otherwise all reads will be used for assembly.

The Virontus pipeline is built using [Nextflow], a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with [Docker] and [Singularity] containers making installation trivial and results highly reproducible.

### Installation

You will need to install [Nextflow] in order to run the Virontus pipeline. 

[Singularity] is recommended for portable and reproducible execution of the pipeline with the `-profile singularity` command-line argument.

#### 1) Install [Nextflow]

If you have [Conda] installed, you can install [Nextflow] with the following command:

```bash
conda install -c bioconda -c conda-forge nextflow
```

#### 2) Install [Singularity][]

Installing [Singularity] is optional but recommended for portability and reproducibility of results. 


#### 3) Install Virontus

Nextflow will automatically download the latest version of Virontus. You can show the Virontus help message with usage information with:

```bash
nextflow run peterk87/nf-virontus --help
```


### Usage

Show usage information with

```bash
nextflow run peterk87/nf-virontus --help
```

You should see the following

```
N E X T F L O W  ~  version 20.01.0
Launching `main.nf` [awesome_pauling] - revision: 9aeb19496b
WARN: DSL 2 IS AN EXPERIMENTAL FEATURE UNDER DEVELOPMENT -- SYNTAX MAY CHANGE IN FUTURE RELEASE
==================================================================
peterk87/nf-virontus   ~  version 1.0.0
==================================================================

  Git info: null - null [null]

Usage:
Given some barcoded and demultiplexed reads, the typical command for running the pipeline is as follows:

  nextflow run peterk87/nf-virontus \
    --reads "reads/*.fastq" \
    --outdir results \
    --ref_fasta refs.fa \
    -profile singularity # recommended to run with Singularity

The above assumes that you have a Centrifuge DB and Kraken2 DB located at
/opt/DB/centrifuge/nt-2018-03-03/nt and /opt/DB/kraken2/standard2,
respectively, OR that you have set $CENTRIFUGE_DB and $KRAKEN2_DB env
variables. It also assumes that you have Singularity installed on your
local machine and will automatically pull and use the Singularity image for
this workflow from Singularity-Hub.org.

NOTE: For best results, please ensure you have Singularity installed prior to running this workflow.(https://sylabs.io/guides/3.3/user-guide/quick_start.html#quick-installation-steps)

Note:
The argument supplied to "--reads" must be quoted if using "*" and other
characters and symbols that could be shell expanded!

Mandatory Options:
  --reads   Input reads directory and pattern (default: "reads/*.fastq")
  --ref_fasta      Reference genomes multiFASTA file (one or more references
                   in a single file) (default: "./refs.fasta")
Amplicon Sequencing Options:
  --bedfile        BED format file with amplicon sequencing primers info (optional).
                   Produced as output from PrimalScheme.
Consensus Generation Options:
  --low_coverage   Low coverage threshold (default=3).
                   Replace consensus sequence positions below this depth
                   threshold with a low coverage character
                   (see --low_cov_char)
  --no_coverage    No coverage threshold (default=0).
                   Replace consensus sequence positions with less than or
                   equal this depth with a no coverage character
                   (see --no_cov_char)
  --low_cov_char   Low coverage character (default="N")
  --no_cov_char    No coverage character (default="-")

Cluster Options:
  --slurm_queue     Name of SLURM queue to run workflow on; use with -profile slurm


Taxonomic Classification Options:
  --centrifuge_db   Path to Centrifuge DB and prefix. If not specified, will
                    try to get from $CENTRIFUGE_DB env variable or see if
                    "/opt/DB/centrifuge/nt-2018-03-03/nt" exists.
                    (default: null)
  --kraken2_db      Path to Kraken2 DB directory. . If not specified, will
                    try to get from $KRAKEN2_DB env variable or see if
                    "/opt/DB/kraken2/standard2" exists.
                    (default: null)
  --taxids          Taxonomic IDs to filter reads by. Multiple taxids should
                    be delimited by commas (`--taxids 1,2,3`). To disable
                    filtering of reads based on taxids, do not provide a
                    value for the `--taxids` argument:
                    `nextflow run ... --taxids --reads ...`
                    (default: 10239 (Viruses))
  --exclude_unclassified_reads  Exclude unclassified reads from taxonomic
                                classification filtered reads (default: false)

De Novo Assembly Options:
  --do_unicycler_assembly       Assemble filtered reads using Unicycler? (default: false)

Other Options:
  --outdir          The output directory where the results will be saved
                    (default: results)
  -w/--work-dir     The temporary directory where intermediate data will be
                    saved (default: ./work)
  -profile          Configuration profile to use. [standard, singularity,
                    conda, slurm] (default 'standard')
  --tracedir        Pipeline run info output directory (default:
                    results/pipeline_info)

Note:
It is recommended that this workflow be executed with Singularity using the
Singularity profile (`-profile singularity`) for maximum reproducibility and
ease of execution on different platforms.
```

#### Preparing your data

It is assumed that your data has been basecalled using the latest version of ONT Guppy (`guppy_basecaller`/`guppy_basecall_server`) and barcode demultiplexed using `guppy_barcoder` with the appropriate settings for the kits used.

After basecalling and demultiplexing, it is recommended that all reads belonging to a particular barcode be concatenated together and optionally renamed to represent the sample to which the reads belong. Virontus will extract the sample name for each input reads FASTQ file from the base filename of the FASTQ file (e.g. sample name will be `sample` from filename `sample1.fastq`).

Below is an example `guppy_barcoder` command for more lenient barcode demultiplexing:

```bash
guppy_barcoder \
  -q 0 \
  --min_score 30 \
  --detect_mid_strand_barcodes \
  --allow_inferior_barcodes \
  --trim_barcodes \
  -i basecalled-reads/ \
  -s demuxed-reads \
  --arrangements_files barcode_arrs_nb12.cfg
```

- `-q 0` to output less files per barcode
- `--min_score 30` for a lower barcode score threshold (default: 60)
- `--detect_mid_strand_barcodes` to detect mid strand barcodes
- `--trim_barcodes` to trim barcodes from read sequences
- `--arrangements_files` to specify the barcodes used

**NOTE:** It's recommended to use the default setting if possible to avoid misassigning reads into the incorrect barcodes.


##### Recommended Steps

1. Basecall reads using Guppy
2. Demultiplex reads using `guppy_barcoder`
3. Concatenate reads belonging to the same barcode into a single file (`cat barcode01/*.fastq > concat-reads/barcode01.fastq`)
4. [Optionally] rename concatenated barcoded reads with appropriate sample name (`mv concat-reads/barcode01.fastq concat-reads/sample1.fastq`)

#### Example

Example command

```bash
$ nextflow run peterk87/nf-virontus \
    -resume \
    -profile singularity \
    --reads "reads/*.fq" \
    --ref_fasta MN908947.3.fa \
    --low_coverage 3 \
    --bedfile nCoV-2019.bed
```

What you will see in the terminal:

```
N E X T F L O W  ~  version 20.01.0
Launching `../main.nf` [ecstatic_davinci] - revision: 9aeb19496b
WARN: DSL 2 IS AN EXPERIMENTAL FEATURE UNDER DEVELOPMENT -- SYNTAX MAY CHANGE IN FUTURE RELEASE
=======================================================
peterk87/nf-virontus v1.0.0
=======================================================
Pipeline Name         : peterk87/nf-virontus
Pipeline Version      : 1.0.0
Run Name              : ecstatic_davinci
Reads                 : reads/*.fq
Ref Sequences FASTA   : MN908947.3.fa
Primer Scheme         : nCoV-2019.bed
Consensus No Coverage : <=0X positions replaced with '-'
Consensus Low Coverage: <3X positions replaced with 'N'
Centrifuge DB         : null
Kraken2 DB            : null
Taxids                : Filtering for taxids belonging to 10239
Unicycler Assembly?   : No
Max Memory            : 256 GB
Max CPUs              : 48
Max Time              : 10d
Output dir            : results
Working dir           : ./work
Container Engine      : singularity
Container             : virontus.simg
Current home          : /home/pkruczkiewicz
Current user          : pkruczkiewicz
Current path          : ./
Script dir            : ./nf-virontus
Config Profile        : standard
Command-Line          : nextflow run peterk87/nf-virontus -profile singularity -resume --reads 'reads/*.fq' --ref_fasta MN908947.3.fa --low_coverage 3 --bedfile nCoV-2019.bed
Nextflow version      : 20.01.0
=========================================
executor >  local (18)
[0a/142458] process > REC2FASTA  [100%] 1 of 1, cached: 1 ✔
[a3/3168c5] process > MAP        [100%] 3 of 3, cached: 3 ✔
[0d/8a698f] process > IVAR_TRIM  [100%] 3 of 3 ✔
[76/f82320] process > MAP_STATS  [100%] 3 of 3 ✔
[cc/de6b36] process > MEDAKA     [100%] 3 of 3 ✔
[74/058b57] process > LONGSHOT   [100%] 3 of 3 ✔
[b4/5ed366] process > BCF_FILTER [100%] 3 of 3 ✔
[a3/ae8e3a] process > CONSENSUS  [ 67%] 2 of 3

Pipeline execution summary
Completed at: 30-Apr-2020 14:00:11
Duration    : 1m 40s
CPU hours   : 0.1 (58.9% cached)
Succeeded   : 18
Cached      : 4
```

Example output file tree structure:


```
results/
├── consensus
│   ├── NB02-MN908947.3.consensus.fasta
│   ├── NB04-MN908947.3.consensus.fasta
│   └── unclassified-MN908947.3.consensus.fasta
├── mapping
│   ├── NB02
│   │   ├── bamfiles
│   │   │   ├── NB02-MN908947.3.bam 
│   │   │   └── NB02-MN908947.3.trim.bam 
│   │   ├── NB02-MN908947.3-depths.tsv
│   │   ├── NB02-MN908947.3.flagstat
│   │   └── NB02-MN908947.3.idxstats
│   ├── NB04
│   │   ├── bamfiles
│   │   │   ├── NB04-MN908947.3.bam 
│   │   │   └── NB04-MN908947.3.trim.bam 
│   │   ├── NB04-MN908947.3-depths.tsv
│   │   ├── NB04-MN908947.3.flagstat
│   │   └── NB04-MN908947.3.idxstats
│   └── unclassified
│       ├── bamfiles
│       │   ├── unclassified-MN908947.3.bam 
│       │   └── unclassified-MN908947.3.trim.bam 
│       ├── unclassified-MN908947.3-depths.tsv
│       ├── unclassified-MN908947.3.flagstat
│       └── unclassified-MN908947.3.idxstats
├── pipeline_info
│   ├── execution_dag.dot
│   ├── execution_report.html
│   ├── execution_timeline.html
│   └── execution_trace.txt
├── refs
│   └── MN908947.3.fa
└── vcf
    ├── NB02-MN908947.3.longshot.filt.vcf
    ├── NB02-MN908947.3.longshot.vcf
    ├── NB02-MN908947.3.medaka.vcf
    ├── NB04-MN908947.3.longshot.filt.vcf
    ├── NB04-MN908947.3.longshot.vcf
    ├── NB04-MN908947.3.medaka.vcf
    ├── unclassified-MN908947.3.longshot.filt.vcf
    ├── unclassified-MN908947.3.longshot.vcf
    └── unclassified-MN908947.3.medaka.vcf
```

### Credits
peterk87/nf-virontus was originally written by Peter Kruczkiewicz.

Bootstrapped with [nf-core/tools](https://github.com/nf-core/tools) `nf-core create`. 

Thank you to the [nf-core/tools](https://github.com/nf-core/tools) team for a great tool for bootstrapping creation of a production ready Nextflow workflows.




[Minimap2]: https://github.com/lh3/minimap2
[Medaka]: https://github.com/nanoporetech/medaka
[Longshot]: https://www.nature.com/articles/s41467-019-12493-y
[Nextflow]: https://www.nextflow.io
[Singularity]: https://sylabs.io/guides/3.5/user-guide/
[Conda]: https://conda.io/
[iVar]: https://github.com/andersen-lab/ivar
[Kraken2]: https://ccb.jhu.edu/software/kraken2/
[Centrifuge]: https://ccb.jhu.edu/software/centrifuge/manual.shtml
[Unicycler]: https://github.com/rrwick/Unicycler
[Docker]: https://www.docker.com/
