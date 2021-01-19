# peterk87/nf-ionampliseq: Usage

## Introduction

This pipeline performs read mapping and variant calling with Thermo Fisher developed open-source tools, [TMAP] and [TVC] to produce more accurate read mapping and variant calling results from Ion Torrent AmpliSeq sequence data. The pipeline also generates a consensus sequence and comprehensive QC stats and results report using MultiQC.

This workflow currently includes several built-in analysis packages for Ion Torrent AmpliSeq sequence data of [CSFV] and [FMDV]. Users can also specify their own AmpliSeq panels, however, these files (reference sequences FASTA and detailed BED file) must be compatible with the Ion Torrent Software Suite including [tmap] and [tvc].

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run peterk87/nf-ionampliseq --input '/path/to/ion-torrent/*.bam' -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

## Help and usage info

The help and usage information for this pipeline can be displayed in your terminal with:

```bash
nextflow run peterk87/nf-ionampliseq --help
```

Which should show a help and usage message like this:

```text
==================================================================
peterk87/nf-ionampliseq   ~  version 1.0.0dev
==================================================================

  Git info: XXX - YYY [ZZZ]

Usage:

The typical command for running the pipeline is as follows:

$ nextflow run peterk87/nf-ionampliseq \
    --input '/path/to/iontorrent/*.bam' \
    --outdir ./results \
    -profile docker # Recommended to run workflow with either Docker or Singularity enabled

Input Options:
  --input           Path to BAM files (e.g. 'ion-torrent/*.bam'). Sample names and AmpliSeq panel will be inferred from the BAM file headers. [Recommended run mode]
  --rundir          Path to Ion Torrent sequencing run containing 'IonCode_*_rawlib.bam' and 'ion_params_00.json' output files.
  --sample_sheet    Sample sheet CSV, TSV, ODS or XLSX file.
  --panel           AmpliSeq panel to run. Choice of 'fmd' or 'csf'. Only needs to be specified with '--sample_sheet'.
  --ref_fasta       Custom AmpliSeq panel reference sequences FASTA.
  --bed_file        Custom AmpliSeq detailed BED file accompanying '--ref_fasta'.

Mash Screen Options:
  --mash_k          Mash sketch kmer size (default: 19)
  --mash_s          Mash sketch number of sketch hashes (default: 10000)

TVC Options:
  --tvc_error_motifs_dir        Directory with Ion Torrent TVC error motifs (default: /home/pkruczkiewicz/sandbox/2020-09-21-fmdv-ion-torrent-weird-gap/ionampliseq/data/tvc-sse)
  --tvc_read_limit              TVC read limit (default: 4000000)
  --tvc_downsample_to_coverage  TVC downsample to at most X coverage (default: X=8000)
  --tvc_min_mapping_qv          TVC min mapping quality value (default: 0)
  --tvc_read_snp_limit          TVC: do not use reads with number of SNPs about this (default: 20)

Cluster Options:
  --slurm_queue       Name of SLURM queue to run workflow on. Must be specified with -profile slurm.
  --slurm_queue_size  Maximum number of Slurm jobs to queue (default: 100)

Other Options:
  --outdir          The output directory where the results will be saved
                    (default: ./results)
  -w/--work-dir     The temporary directory where intermediate data will be
                    saved (default: /home/pkruczkiewicz/sandbox/2020-09-21-fmdv-ion-torrent-weird-gap/ionampliseq/test/derp/work)
  -profile          Configuration profile to use. [standard, singularity,
                    conda, slurm] (default 'standard')
  --tracedir        Pipeline run info output directory (default:
                    ./results/pipeline_info)
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull peterk87/nf-ionampliseq
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [peterk87/nf-ionampliseq releases page](https://github.com/peterk87/nf-ionampliseq/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](https://docker.com/)
  * Pulls software from Docker Hub: [`peterk87/nf-ionampliseq`](https://hub.docker.com/r/peterk87/nf-ionampliseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
  * Pulls software from Docker Hub: [`peterk87/nf-ionampliseq`](https://hub.docker.com/r/peterk87/nf-ionampliseq/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

#### Custom resource requests

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

Whilst these default requirements will hopefully work for most people with most data, you may find that you want to customise the compute resources that the pipeline requests. You can do this by creating a custom config file. For example, to give the workflow process `TVC` 32GB of memory, you could use the following config:

```nextflow
process {
  withName: TVC {
    memory = 32.GB
  }
}
```

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

<!-- TODO: uncomment when pipeline added to nf-core

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).
-->

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

#### Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
