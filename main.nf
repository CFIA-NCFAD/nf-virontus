#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { checkFileExists; checkKraken2Db; check_sample_sheet } from './lib/helpers'

def json_schema = "$projectDir/nextflow_schema.json"
// Show help message if --help specified
if (params.help){
  def command = "nextflow run CFIA-NCFAD/nf-virontus --input samplesheet.csv --genome 'MN908947.3' --artic_v4 -profile docker"
  log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
  exit 0
}

//=============================================================================
// CHECK PARAMS
//=============================================================================

if (params.validate_params) {
  NfcoreSchema.validateParameters(params, json_schema, log)
}

if (workflow.profile == 'slurm' && params.slurm_queue == "") {
  log.error "You must specify a valid SLURM queue (e.g. '--slurm_queue <queue name>' (see `\$ sinfo` output for available queues)) to run this workflow with the 'slurm' profile!"
  exit 1
}

// if (params.kraken2_db) {
//   checkKraken2Db(params.kraken2_db)
// }

genome = params.genome
fasta = params.fasta
gff = params.gff
primer_bed = params.primer_bed
// Viral reference files
if (params.scov2) {
  genome = 'MN908947.3'
}
if (params.genomes && genome && !params.genomes.containsKey(genome)) {
  exit 1, "The provided genome '${genome}' is not available in the Genome file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}
if (!fasta || fasta == null) {
  fasta = genome ? params.genomes[genome].fasta : false
} else {
  log.info "Using user-specified reference FASTA '${fasta}'"
}
if (!gff || gff == null) { 
  gff = genome ? params.genomes[genome].gff : false
} else {
  log.info "Using user-specified reference GFF '${gff}'"
}
if (fasta) {
  file(fasta, checkIfExists: true)

  lastPath = fasta.lastIndexOf(File.separator)
  lastExt = fasta.lastIndexOf(".")
  fasta_base = fasta.substring(lastPath+1)
  index_base = fasta.substring(lastPath+1,lastExt)
  if (fasta.endsWith('.gz')) {
      fasta_base = fasta.substring(lastPath+1,lastExt)
      index_base = fasta_base.substring(0,fasta_base.lastIndexOf("."))
  }
} else {
  exit 1, "Viral genome fasta file not specified!"
}

if (genome == 'MN908947.3' && !primer_bed) {
  if (params.artic_v3) {
    primer_bed = params.genomes[genome]['primer_schemes']['artic_v3']
  } else if (params.artic_v4) {
    primer_bed = params.genomes[genome]['primer_schemes']['artic_v4']
  } else if (params.artic_v4_1) {
    primer_bed = params.genomes[genome]['primer_schemes']['artic_v4_1']
  } else if (params.freed) {
    primer_bed = params.genomes[genome]['primer_schemes']['freed']
  } else if (params.neb_primers) {
    primer_bed = params.genomes[genome]['primer_schemes']['neb_primers']
  } else {
    log.warn "Using SARS-CoV-2 ${genome} as reference genome with no primer BED file. Primers will not be trimmed since they are not specified. You can specify a custom primer scheme in BED file format or use one of the built-in primer schemes: ${params.genomes[genome].primer_schemes.keySet().join(', ')}"
  }
}

if (genome == 'MN908947.3') {
  nextclade_dataset = 'sars-cov-2'
} else {
  nextclade_dataset = params.nextclade_dataset
}

if (params.clair3_user_variant_model) {
  ch_user_clair3_model = file(params.clair3_user_variant_model, checkIfExists: true)
}

//=============================================================================
// LOG PARAMS SUMMARY
//=============================================================================

def summary_params = NfcoreSchema.params_summary_map(workflow, params, json_schema)
log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)


//=============================================================================
// PROCESSES
//=============================================================================
include { CHECK_SAMPLE_SHEET } from './modules/local/check_sample_sheet'
include { CAT_FASTQ } from './modules/local/cat_fastq'
include { MINIMAP2 } from './modules/local/minimap2'
include { PRIMER_BED_TO_AMPLICON_BED } from './modules/local/primer_bed_to_amplicon_bed'
include { MOSDEPTH } from './modules/local/mosdepth'
include { IVAR_TRIM } from './modules/local/ivar'
include { CLAIR3 } from './modules/local/clair3'

// include { MULTIQC_TSV_FROM_LIST as READ_COUNT_FAIL_TSV        } from '../modules/local/multiqc_tsv_from_list'
// include { MULTIQC_TSV_FROM_LIST as READ_COUNT_PASS_TSV        } from '../modules/local/multiqc_tsv_from_list'

include { BCFTOOLS_STATS as BCFTOOLS_STATS_PRE_FILTER } from './modules/local/bcftools_stats'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_POST_FILTER } from './modules/local/bcftools_stats'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_MINOR } from './modules/local/bcftools_filter'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_MAJOR } from './modules/local/bcftools_filter'
include { VCF_FILTER_FRAMESHIFT as VCF_FILTER_FRAMESHIFT_MINOR } from './modules/local/vcf_filter_frameshift'
include { VCF_FILTER_FRAMESHIFT as VCF_FILTER_FRAMESHIFT_MAJOR } from './modules/local/vcf_filter_frameshift'

include { SNPEFF_BUILD } from './modules/local/snpeff_build'
include { SNPEFF_ANN } from './modules/local/snpeff_ann'
include { TABIX as TABIX_SNPEFF } from './modules/local/tabix'
include { TABIX as TABIX_CONSENSUS } from './modules/local/tabix'
include { SNPSIFT } from './modules/local/snpsift'
include { SIMPLER_SNPSIFT } from './modules/local/simpler_snpsift'
include { BCFTOOLS_CONSENSUS } from './modules/local/bcftools_consensus'
include { COVERAGE_PLOT                      } from './modules/local/plots'
include { MULTIQC; CONSENSUS_TO_MULTIQC_HTML } from './modules/local/multiqc'

include { PYCOQC } from './modules/local/pycoqc'
include { NANOPLOT } from './modules/local/nanoplot'
// include { KRAKEN2_RUN } from './modules/local/kraken2'

include { PREPARE_FASTA_FOR_PANGOLIN; PANGOLIN; PANGOLIN_SUMMARY_FOR_MULTIQC } from './modules/local/pangolin'
include { NEXTCLADE_DATASETGET; NEXTCLADE_RUN } from './modules/local/nextclade'
include { KRAKEN2_PREPAREINDEX } from './modules/local/kraken2'

include { FGBIO_CLIPBAM } from './modules/local/fgbio_clipbam'
include { SAMTOOLS_FASTQ } from './modules/local/samtools_fastq'
include { SAMTOOLS_SORT_FGBIO } from './modules/local/samtools_sort_fgbio'

//=============================================================================
// MAIN WORKFLOW
//=============================================================================
workflow {

  ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)
  ch_versions = Channel.empty()

  if (params.kraken2_db) {
    ch_kraken2_index = file(params.kraken2_db)
    KRAKEN2_PREPAREINDEX(ch_kraken2_index)
    KRAKEN2_PREPAREINDEX.out.view()
  }

  // If sample sheet table provided, 
  //   - validate and write to CSV sample sheet
  //   - use 'check_sample_sheet' method to fetch files for each sample 
  //     (e.g. all FASTQ files within a Guppy barcoding output directory like 
  //     'barcode01/')
  //   - concatenate multiple FASTQs and gzip compress into single fastq.gz 
  //     per sample
  CHECK_SAMPLE_SHEET(Channel.from(file(params.input, checkIfExists: true)))
  CHECK_SAMPLE_SHEET.out
    .splitCsv(header: ['sample', 'reads'], sep: ',', skip: 1)
        // "reads" can be path to file or directory
    .map { [it.sample, it.reads] }
    // group by sample name to later merge all reads for that sample
    .groupTuple(by: 0)
    // collect all uncompressed and compressed FASTQ reads into 2 lists
    // and count number of reads for sample
    .map { sample, reads ->
      // uncompressed FASTQ list
      def fq = []
      // compressed FASTQ list
      def fqgz = []
      // read count
      def count = 0
      for (f in reads) {
        f = file(f)
        if (f.isFile() && f.getName() ==~ /.*\.(fastq|fq)(\.gz)?/) {
          if (f.getName() ==~ /.*\.gz/) {
            fqgz << f
          } else {
            fq << f
          }
          continue
        }
        if (f.isDirectory()) {
          // only look for FQ reads in first level of directory
          for (x in f.listFiles()) {
            if (x.isFile() && x.getName() ==~ /.*\.(fastq|fq)(\.gz)?/) {
              if (x.getName() ==~ /.*\.gz/) {
                fqgz << x
              } else {
                fq << x
              }
            }
          }
        }
      }
      for (x in fq) {
        count += x.countFastq()
      }
      for (x in fqgz) {
        count += x.countFastq()
      }
      return [ sample, fqgz, fq, count ]
    }
    .set { ch_input_sorted }

  // ch_input_sorted
  //   .branch { sample, fqgz, fq, count  ->
  //     pass: count >= params.min_sample_reads
  //       pass_sample_reads[sample] = count
  //       return [ "$sample\t$count" ]
  //     fail: count < params.min_sample_reads
  //       fail_sample_reads[sample] = count
  //       return [ "$sample\t$count" ]
  //   }
  //   .set { ch_pass_fail_read_count }

  // // Report samples which have reads count < min_sample_reads
  // READ_COUNT_FAIL_TSV(
  //   ch_pass_fail_read_count.fail.collect(),
  //   ['Sample', 'Read count'],
  //   'fail_read_count_samples'
  // )
  // // Report samples which have reads count >= min_sample_reads
  // READ_COUNT_PASS_TSV(
  //   ch_pass_fail_read_count.pass.collect(),
  //   ['Sample', 'Read count'],
  //   'pass_read_count_samples'
  // )

  // Keep samples which have reads count  > min_sample_reads for downstream analysis
  // Re-arrange channels to have meta map of information for sample
  ch_input_sorted
    .filter { it[-1] >= params.min_sample_reads }
    .map { sample, fqgz, fq, count -> [ [id: sample], fqgz, fq ] }
    .set { ch_reads }

  CAT_FASTQ(ch_reads)

  CAT_FASTQ.out.reads.set {ch_reads_cat}

  if (primer_bed) {
    ch_bed = Channel.value(file(primer_bed))
  } else {
    ch_bed = ch_dummy_file
  }

  if (params.sequencing_summary && !params.skip_pycoqc) {
    PYCOQC(ch_sequencing_summary)
    ch_versions = ch_versions.mix(PYCOQC.out.versions)
  }

  // Nanoplot read level QC
  if (!params.skip_nanoplot) {
    NANOPLOT(ch_reads_cat)
    ch_versions = ch_versions.mix(NANOPLOT.out.versions)
  }

  // reference fasta into channel
  ch_fasta = Channel.value(file(fasta))
  // if reference GFF specified, create SnpEff DB
  if (gff) {
    SNPEFF_BUILD(file(fasta), file(gff))
    ch_versions = ch_versions.mix(SNPEFF_BUILD.out.versions.first().ifEmpty(null))
  }

  // Optional host subtraction with Kraken2
  // if (params.kraken2_db && !params.skip_kraken2) {
  //   KRAKEN2_RUN(file(params.kraken2_db), ch_reads_cat)
  //   if (params.subtract_host) {
  //     ch_reads_cat = KRAKEN2_RUN.out.unclassified
  //   }
  // }
  // Map reads to reference
  MINIMAP2(ch_reads_cat.combine(ch_fasta))
  ch_versions = ch_versions.mix(MINIMAP2.out.versions.first().ifEmpty(null))
  // Trim primer sequences from read alignments if primer scheme BED file provided 
  if (primer_bed) {
    IVAR_TRIM(MINIMAP2.out.bam, ch_bed)
    ch_versions = ch_versions.mix(IVAR_TRIM.out.versions.first().ifEmpty(null))
    ch_bam = IVAR_TRIM.out.bam
    ch_depths = IVAR_TRIM.out.depth
  } else {
    ch_bam = MINIMAP2.out.bam
    ch_depths = MINIMAP2.out.depth
  }

  if (params.output_hardclipped_reads) {
    SAMTOOLS_SORT_FGBIO(ch_bam)
    FGBIO_CLIPBAM(SAMTOOLS_SORT_FGBIO.out.bam, ch_fasta)
    ch_versions = ch_versions.mix(FGBIO_CLIPBAM.out.versions.first().ifEmpty(null))

    SAMTOOLS_FASTQ(FGBIO_CLIPBAM.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions.first().ifEmpty(null))
  }

  PRIMER_BED_TO_AMPLICON_BED(ch_bed)
  MOSDEPTH(ch_bam, PRIMER_BED_TO_AMPLICON_BED.out, params.mosdepth_window_size)
  ch_versions = ch_versions.mix(MOSDEPTH.out.versions.ifEmpty(null))

  if (params.clair3_user_variant_model) {
    CLAIR3(ch_bam, ch_fasta, ch_user_clair3_model)
  } else {
    CLAIR3(ch_bam, ch_fasta, [])
  }
  ch_versions = ch_versions.mix(CLAIR3.out.versions.first().ifEmpty(null))

  BCFTOOLS_FILTER_MINOR(CLAIR3.out.vcf, params.minor_allele_fraction)
  ch_versions = ch_versions.mix(BCFTOOLS_FILTER_MINOR.out.versions.first().ifEmpty(null))
  BCFTOOLS_FILTER_MAJOR(CLAIR3.out.vcf, params.major_allele_fraction)
  VCF_FILTER_FRAMESHIFT_MINOR(BCFTOOLS_FILTER_MINOR.out.vcf)
  VCF_FILTER_FRAMESHIFT_MAJOR(BCFTOOLS_FILTER_MAJOR.out.vcf)
  BCFTOOLS_STATS_PRE_FILTER(VCF_FILTER_FRAMESHIFT_MINOR.out, ch_fasta)
  ch_versions = ch_versions.mix(BCFTOOLS_STATS_PRE_FILTER.out.versions.first().ifEmpty(null))
  BCFTOOLS_STATS_POST_FILTER(VCF_FILTER_FRAMESHIFT_MAJOR.out, ch_fasta)

  if (gff) {
    SNPEFF_ANN(
      VCF_FILTER_FRAMESHIFT_MINOR.out,
      SNPEFF_BUILD.out.db,
      SNPEFF_BUILD.out.config,
      ch_fasta
    )
    ch_versions = ch_versions.mix(SNPEFF_ANN.out.versions.first().ifEmpty(null))

    TABIX_SNPEFF(SNPEFF_ANN.out.vcf)
    ch_versions = ch_versions.mix(TABIX_SNPEFF.out.versions.first().ifEmpty(null))

    SNPSIFT(TABIX_SNPEFF.out.tbi)
    ch_versions = ch_versions.mix(SNPSIFT.out.versions.first().ifEmpty(null))
    
    SIMPLER_SNPSIFT(SNPSIFT.out.txt)
    ch_snpeff = SNPEFF_ANN.out.csv | collect
  } else {
    ch_snpeff = Channel.from([])
  }
  TABIX_CONSENSUS(VCF_FILTER_FRAMESHIFT_MAJOR.out)
  BCFTOOLS_CONSENSUS(
    TABIX_CONSENSUS.out.tbi.join(MOSDEPTH.out.per_base_bed),
    ch_fasta,
    params.low_coverage
  )
  ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions.first().ifEmpty(null))
  if (!params.skip_coverage_plot) {
    COVERAGE_PLOT(
      VCF_FILTER_FRAMESHIFT_MAJOR.out.join(ch_depths),
      ch_fasta
    )
    // ch_versions = ch_versions.mix(COVERAGE_PLOT.out.versions.first().ifEmpty(null))
  }
  if (genome == 'MN908947.3') {
    if (params.tree_extra_fasta) {
      ch_extra_fasta_with_sample_name = Channel.fromPath(params.tree_extra_fasta).splitFasta(file: true) \
        | map {
          f = file(it)
          m = f.text =~ /^>(\S+).*/
          sample = m[0][1]
          sample = sample.replaceAll(/[^\w\-]/, "_")
          [sample, it]
        }
      ch_fasta_for_pangolin = BCFTOOLS_CONSENSUS.out.consensus | mix(ch_extra_fasta_with_sample_name)
    } else {
      ch_fasta_for_pangolin = BCFTOOLS_CONSENSUS.out.consensus
    }
    PREPARE_FASTA_FOR_PANGOLIN(ch_fasta_for_pangolin)
    PANGOLIN(
      PREPARE_FASTA_FOR_PANGOLIN.out.collectFile(name: "sequences.fasta", newLine: true)
    )
    ch_versions = ch_versions.mix(PANGOLIN.out.versions.first().ifEmpty(null))
    PANGOLIN_SUMMARY_FOR_MULTIQC(PANGOLIN.out.report)
    ch_pangolin = PANGOLIN.out.report
    ch_pangolin_mqc = PANGOLIN_SUMMARY_FOR_MULTIQC.out

    NEXTCLADE_DATASETGET(params.nextclade_dataset)
    ch_fastas = BCFTOOLS_CONSENSUS.out.consensus.map { it[1] }.collect()
    NEXTCLADE_RUN(ch_fastas, NEXTCLADE_DATASETGET.out.nextclade_dataset)
    ch_versions = ch_versions.mix(NEXTCLADE_RUN.out.versions.first().ifEmpty(null))
  } else {
    ch_pangolin = Channel.from([])
    ch_pangolin_mqc = Channel.from([])
  }
}

//=============================================================================
// INTROSPECTION
//=============================================================================
workflow.onComplete {
    println """
    Pipeline execution summary
    ---------------------------
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Results Dir  : ${file(params.outdir)}
    Work Dir     : ${workflow.workDir}
    Exit status  : ${workflow.exitStatus}
    Error report : ${workflow.errorReport ?: '-'}
    """.stripIndent()
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
