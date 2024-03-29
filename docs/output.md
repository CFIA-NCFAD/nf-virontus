# CFIA-NCFAD/nf-virontus: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [Concatenated FASTQ Reads](#concatenated-fastq-reads) - Concatenation and Gzip compression of FASTQ reads
* [Read Mapping](#read-mapping) - Read mapping using [Minimap2][] and read mapping stats calculation with [Samtools][]
* [Mosdepth](#mosdepth) - Coverage stats calculated by [Mosdepth][]
* [Variant Calling](#variant-calling) - Variant calling using [Clair3] followed by filtering for major/minor variants is performed for SnpEff variant effect analysis, consensus sequence generation and high confidence variant statistics for MultiQC report.
* [SnpEff](#snpeff) - SnpEff variant effect analysis on major and minor variants.
* [Consensus Sequence](#consensus-sequence) - Consensus sequence with `N` masking of low/no coverage positions.
* [Pangolin](#pangolin) - SARS-CoV-2 global lineage assignment using [Pangolin][] if using SARS-CoV-2 Wuhan-Hu-1 MN908947.3 as reference.
* [Coverage Plots](#coverage-plots) - Coverage plots with/without low/no coverage and/or variants highlighted with linear and log10 scaling of y-axis depth values.
* [MultiQC](#multiqc) - Aggregate report describing results from the whole pipeline
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Concatenated FASTQ Reads

When a user supplies a sample sheet CSV and read data paths point to directories (e.g. `fastq_pass/barcode01`), then the FASTQ files are concatenated together, Gzip compressed with [pigz] and output with the filename `${sample}.fastq.gz`.

**Output files:**

* `reads/`
  * `*.fastq.gz`: Gzip compressed concatenated FASTQ format read sequences extracted for each sample if `--save_cat_reads` specified when running the workflow.

## Read Mapping

[Minimap2][] is used for read mapping the Nanopore reads against the reference genome.

[Samtools][] is used to compute read mapping statistics and metrics including coverage depth at each position of the reference genome.

**Output files:**

* `mapping/${sample}`
  * `*.{fa,fasta}`: Reference genome FASTA.
  * `*.bam`: BAM file with read alignment of sample reads against reference.
  * `*.bam.bai`: BAM file index (`samtools index` output).
  * `*.flagstat`: The number of alignments for each FLAG type.
  * `*.idxstats`: Alignment summary statistics.
  * `*.stats`: Comprehensive read alignment statistics.
* `samtools/depth/`
  * `*-depths.tsv`: Tab-delimited table of read mapping coverage depth at each position in the reference containing 3 columns: reference genome ID, position and coverage depth.

> **NB:** Stats and metrics computed by Samtools are used to generate multiple fields in the MultiQC general statistics table and to generate multiple plots in the MultiQC report.

## Mosdepth

[Mosdepth][] is used to compute coverage depth and breadth statistics. Mosdepth results for each sample are visualized and summarized in the MultiQC report.

**Output files:**

* `mosdepth/`
  * `*.mosdepth.global.dist.txt`: contains, a cumulative distribution indicating the proportion of total bases  that were covered for at least a given coverage value. It does this for each chromosome, and for the whole genome.
  * `*.mosdepth.region.dist.txt`: contains, a cumulative distribution indicating the proportion of 200 base window regions that were covered for at least a given coverage value. It does this for each chromosome, and for the whole genome.
  * `*.mosdepth.summary.txt`: summary read alignment depth statistics (reference sequence length, total bases, mean/min/max).
  * `*.per-base.bed.gz`: 4 column BED file with per base coverage. Columns in order contain: reference ID; 0-base start index; 0-base end index; coverage depth.
  * `*.regions.bed.gz`: 4 column BED file with 200 base window regions mean coverage. Columns in order contain: reference ID; 0-base start index; 0-base end index; coverage depth.
  * `*.gz.csi`: Table index files for rapid random read access.

## Variant Calling

[Clair3][] performs variant calling using deep learning neural network models trained on a variety of ONT output.

Variants detected by [Clair3] are filtered with [Bcftools][] by allele fraction for the subset that qualify as minor or major variants (default: 0.25 and 0.75 allele fraction). Bcftools is also used to report the number of SNPs, MNPs and indels detected as well as transitions/transversions (Ts/Tv). Bcftools statistics are reported in the MultiQC report general stats table.

By default, the minimum allele fraction (AF) or frequency an alternate allele must be observed is 0.25, otherwise the allele is discarded from downstream analysis.

By default, the maximum allele fraction or frequency an alternate allele must be observed to be incorporated into the consensus sequence is 0.75.

Indel variants that potentially introduce frameshift mutations are also filtered out by default.

**Output files:**

* `variants/`
  * `${sample}.clair3/`: Clair3 output directory
  * `${sample}.0.75AF.filt.no_fs.vcf`: 75% AF filtered VCF with potential frameshift mutations removed
  * `${sample}.0.25AF.filt.no_fs.vcf`: 25% AF filtered VCF with potential frameshift mutations removed
  * `${sample}.0.75AF.filt.vcf`: 75% AF filtered VCF
  * `${sample}.0.25AF.filt.vcf`: 25% AF filtered VCF
  * `${sample}.clair3.fix.vcf.gz`: unfiltered Clair3 VCF

## SnpEff

[SnpEff][] is used to annotate variants with variant effect and impact.

[SnpSift][] is used to filter variants by variant effect and impact and to convert the VCF file to a tab-delimited table format.

**Output files:**

* `variants/snpeff/`
  * `*.snpeff.vcf`: SnpEff annotated VCF file
  * `*.snpeff.vcf.gz`: Gzip compressed SnpEff annotated VCF file
  * `*.snpeff.vcf.gz.tbi`: Tabix index file for rapid random read access
  * `*.snpsift.vcf`: SnpSift annotated VCF file
  * `*.snpsift.vcf.gz`: Gzip compressed SnpSift annotated VCF file
  * `*.snpsift.vcf.gz.tbi`: Tabix index file for rapid random read access
  * `*.snpsift.tsv`: SnpSift annotated VCF file converted to tab-delimited table format
  * `*.snpsift.tsv.gz`: Gzip compressed SnpSift annotated VCF file converted to tab-delimited table format
  * `*.snpsift.tsv.gz.tbi`: Tabix index file for rapid random read access

## Consensus Sequence

A consensus sequence is constructed from a coverage depth masked reference sequence (low/no coverage reference sequence positions are replaced with `N`) and selected variants in a VCF file using [Bcftools][].

**Output files:**

* `variants/consensus/`
  * `*.consensus.fasta`: Majority consensus sequence with low/no coverage positions masked with `N` (by default).

## Pangolin

If the workflow is being run with the SARS-CoV-2 reference genome Wuhan-Hu-1 [MN908947.3](https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3/), then [Pangolin][] is used to assign a global SARS-CoV-2 lineage to each sample.

**Output files:**

* `pangolin/`
  * `*.pangolin.lineage_report.csv`: Pangolin lineage assignment report for each sample
  * `logs/`
    * `*.log`: Pangolin logs

> **NB:** In order to run Pangolin with this workflow, you must run the workflow with the `docker` or `singularity` profile (i.e. `nextflow run CFIA-NCFAD/nf-virontus -profile docker/singularity ...`) and you must run the workflow in online mode.

## Coverage Plots

Coverage plots (reference position vs coverage depth) are generated for each sample from the read alignment depths at
each position and filtered variants detected. Several plots are created for each sample with and without low/no coverage
and/or variants highlighted with linear and log10 scaling of y-axis depth values. Coverage plots figures are generated
using a simple Python script that uses the [Matplotlib][] and [seaborn][] plotting libraries.

**Output files:**

* `plots/`
  * `coverage_plot-*.pdf`: coverage plot with low/no coverage regions highlighted in red for no coverage and yellow for low coverage (<= 3X). Linear scaling of y-axis depth.
  * `coverage_plot-*-with_variants.pdf`: coverage plot with low/no coverage regions highlighted in red for no coverage and yellow for low coverage (<= 3X) and TVC variants highlighted. Linear scaling of y-axis depth.
  * `coverage_plot-*-no_low_cov_highlighting.pdf`: coverage plot no highlighting of no/low coverage regions.
  * `coverage_plot-*-log_scale.pdf`: coverage plot with log10 scaling of y-axis depth.

> **NB:** If you're interested in the coverage plot with no annotations, look for the `${sample}-no_low_cov_highlighting.pdf` or `${sample}-no_low_cov_highlighting-log_scale.pdf` if you want log10 scaling of the depth values in the y-axis.

## MultiQC

[MultiQC] is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability.

For more information about how to use MultiQC reports, see [https://multiqc.info](https://multiqc.info).

**Output files:**

* `multiqc/`  
  * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  * `multiqc_plots/`: directory containing static images from the report in various formats.

> **NB:** All consensus sequence FASTA files will be embedded within the MultiQC HTML report and can be downloaded from it.

## Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

**Output files:**

* `pipeline_info/`
  * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
  * Documentation for interpretation of results in HTML format: `results_description.html`.

[Bcftools]: https://samtools.github.io/bcftools/bcftools.html
[Clair3]: https://github.com/HKU-BAL/Clair3
[Matplotlib]: https://matplotlib.org/
[Minimap2]: https://github.com/lh3/minimap2
[Mosdepth]: https://github.com/brentp/mosdepth
[MultiQC]: http://multiqc.info
[Pangolin]: https://github.com/cov-lineages/pangolin/
[pigz]: https://www.zlib.net/pigz/
[Samtools]: https://www.htslib.org/
[seaborn]: https://seaborn.pydata.org/
[SnpEff]: https://pcingola.github.io/SnpEff/
[SnpSift]: https://pcingola.github.io/SnpEff/ss_introduction/
