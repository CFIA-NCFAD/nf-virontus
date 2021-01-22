#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
  helpMessage;
  checkFileExists;
  checkKraken2Db;
  checkTaxids;
  check_sample_sheet;
} from './lib/helpers'

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}

if (workflow.profile == 'slurm' && params.slurm_queue == "") {
  log.error "You must specify a valid SLURM queue (e.g. '--slurm_queue <queue name>' (see `\$ sinfo` output for available queues)) to run this workflow with the 'slurm' profile!"
  exit 1
}

//=============================================================================
// PROCESSES
//=============================================================================

include {
  SOFTWARE_VERSIONS
} from './processes/docs'
include {
  CHECK_SAMPLE_SHEET;
  CAT_FASTQ
} from './processes/misc'
include {
  MAP
} from './processes/mapping'
include {
  MOSDEPTH_GENOME
} from './processes/mosdepth'
include {
  IVAR_TRIM
} from './processes/ivar'
include {
  MEDAKA_LONGSHOT;
  VARIANT_FILTER as VARIANT_FILTER_MAJOR;
  VARIANT_FILTER as VARIANT_FILTER_MINOR;
  BCFTOOLS_STATS as BCFTOOLS_STATS_PRE_FILTER;
  BCFTOOLS_STATS as BCFTOOLS_STATS_POST_FILTER
} from './processes/variants'
include {
  MAKE_SNPEFF_DB;
  SNPEFF
} from './processes/snpeff'
include {
  CONSENSUS
} from './processes/consensus'
include {
  COVERAGE_PLOT
} from './processes/plots'
include {
  MULTIQC;
  CONSENSUS_TO_MULTIQC_HTML
} from './processes/multiqc'


//=============================================================================
// MAIN WORKFLOW
//=============================================================================
workflow {
  if (params.kraken2_db) {
    taxids = checkTaxids(params.taxids)
    checkKraken2Db(params.kraken2_db)
  }

  // Has the run name been specified by the user?
  //  this has the bonus effect of catching both -name and --name
  custom_runName = params.name
  if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
    custom_runName = workflow.runName
  }

  // courtesy of nf-core/viralrecon (https://github.com/nf-core/viralrecon)
  genome = params.genome
  fasta = params.fasta
  gff = params.gff
  bed = params.bed
  // Viral reference files
  if (params.scov2) {
    genome = 'MN908947.3'
    fasta = params.genomes[genome].fasta
    gff = params.genomes[genome].gff
  }
  if (params.genomes && genome && !params.genomes.containsKey(genome)) {
    exit 1, "The provided genome '${genome}' is not available in the Genome file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
  }
  if (!fasta || fasta == null) {
    fasta = genome ? params.genomes[genome].fasta : false
  }
  if (!gff || gff == null) { 
    gff = genome ? params.genomes[genome].gff : false
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

  if (genome == 'MN908947.3' && !bed) {
    if (params.articv3) {
      bed = params.genomes[genome]['schemes']['ARTIC_V3']
    } else if (params.freed) {
      bed = params.genomes[genome]['schemes']['Freed']
    } else {
      log.warn "Using SARS-CoV-2 ${genome} as reference genome with no primer BED file. Primers will not be trimmed since they are not specified. You can specify a custom primer scheme in BED file format or use one of the built-in primer schemes: ${params.genomes[genome]['schemes'].keySet().join(', ')}"
    }
  }

  //=============================================================================
  // LOG EXECUTION START PARAMS
  //=============================================================================

  def summary = [:]
  summary['Pipeline Name']       = workflow.manifest.name
  summary['Pipeline Version']    = workflow.manifest.version
  if (workflow.revision) summary['Pipeline Release'] = workflow.revision
  summary['Run Name']            = custom_runName ?: workflow.runName
  summary['Reads']               = params.reads
  if (genome) summary['Ref Genome'] = genome
  summary['Ref Genome FASTA'] = fasta
  if (gff) summary['Ref Genome GFF'] = gff
  if (bed) summary['Primer Scheme'] = bed
  summary['Major Allele Fraction'] = params.major_allele_fraction
  summary['Minor Allele Fraction'] = params.minor_allele_fraction
  summary['Consensus No Coverage'] = "<=${params.no_coverage}X positions replaced with '${params.no_cov_char}'"
  summary['Consensus Low Coverage'] = "<${params.low_coverage}X positions replaced with '${params.low_cov_char}'"
  if (params.kraken2_db) summary['Kraken2 DB'] = params.kraken2_db
  summary['Max Memory']       = params.max_memory
  summary['Max CPUs']         = params.max_cpus
  summary['Max Time']         = params.max_time
  summary['Current home']     = "$HOME"
  summary['Current user']     = "$USER"
  summary['Current path']     = "$PWD"
  summary['Working dir']      = workflow.workDir
  summary['Output dir']       = params.outdir
  summary['Script dir']       = workflow.projectDir
  summary['Config Profile']   = workflow.profile
  summary['Command-Line']     = workflow.commandLine
  summary['Nextflow version'] = workflow.nextflow.version
  if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
  log.info """=======================================================
    ${workflow.manifest.name} v${workflow.manifest.version}
    =======================================================""".stripIndent()
  log.info summary.collect { k,v -> "${k.padRight(22)}: $v" }.join("\n")
  log.info "========================================="

  if (params.reads) {
    ch_reads = Channel.fromPath(params.reads)
      .map {
        f = file(it)
        filename = f.getName()
        last_ext = filename.lastIndexOf(".")
        index_base = filename.substring(0, last_ext)
        if (filename.endsWith('.gz')) {
          fastq_base = filename.substring(0, last_ext)
          index_base = fastq_base.substring(0, fastq_base.lastIndexOf("."))
        } 
        [index_base, it]
      }
  } else {
    ch_reads = Channel.from([])
  }
  if (params.sample_sheet) {
    // If sample sheet table provided, 
    //   - validate and write to CSV sample sheet
    //   - use 'check_sample_sheet' method to fetch files for each sample 
    //     (e.g. all FASTQ files within a Guppy barcoding output directory like 
    //     'barcode01/')
    //   - concatenate multiple FASTQs and gzip compress into single fastq.gz 
    //     per sample
    //   - add in reads specified by '--reads' CLI parameter
    Channel.from(file(params.sample_sheet, checkIfExists: true)) \
      | CHECK_SAMPLE_SHEET \
      | splitCsv(header: ['sample', 'reads'], sep: ',', skip: 1) \
      | map { check_sample_sheet(it) } \
      | CAT_FASTQ \
      | mix(ch_reads) \
      | set {ch_reads}
  }
  // reference fasta into channel
  ch_fasta = Channel.value(file(fasta))
  // if reference GFF specified, create SnpEff DB
  if (gff) {
    MAKE_SNPEFF_DB(Channel.value([index_base, file(fasta), file(gff)]))
  }
  // Map reads to reference
  ch_reads | combine(ch_fasta) | MAP

  // Trim primer sequences from read alignments if primer scheme BED file provided 
  if (bed) {
    IVAR_TRIM(Channel.value(file(bed)), MAP.out.bam)
    ch_bam = IVAR_TRIM.out.bam
    ch_depths = IVAR_TRIM.out.depths
  } else {
    ch_bam = MAP.out.bam
    ch_depths = MAP.out.depths
  }
  // only interested in sample name [0] and bam [2] for mosdepth
  ch_bam | map {[it[0], it[2]]} | MOSDEPTH_GENOME
  // variant calling and stats
  MEDAKA_LONGSHOT(ch_bam)
  VARIANT_FILTER_MINOR(MEDAKA_LONGSHOT.out, params.minor_allele_fraction) \
    | BCFTOOLS_STATS_PRE_FILTER
  if (gff) {
    SNPEFF(VARIANT_FILTER_MINOR.out, MAKE_SNPEFF_DB.out)
    ch_snpeff = SNPEFF.out.csv | collect
  } else {
    ch_snpeff = Channel.from([])
  }
  VARIANT_FILTER_MAJOR(MEDAKA_LONGSHOT.out, params.major_allele_fraction)
  VARIANT_FILTER_MAJOR.out | BCFTOOLS_STATS_POST_FILTER
  VARIANT_FILTER_MAJOR.out | join(ch_depths) | (CONSENSUS & COVERAGE_PLOT)

  CONSENSUS.out | map { it[1] } | collect | CONSENSUS_TO_MULTIQC_HTML

  if (genome == 'MN908947.3') {
    include {
      PREPARE_FASTA_FOR_PANGOLIN;
      PANGOLIN;
      PANGOLIN_SUMMARY_FOR_MULTIQC
    } from './processes/pangolin'
    if (params.tree_extra_fasta) {
      ch_extra_fasta_with_sample_name = Channel.fromPath(params.tree_extra_fasta).splitFasta(file: true) \
        | map {
          f = file(it)
          m = f.text =~ /^>(\S+).*/
          sample = m[0][1]
          sample = sample.replaceAll(/[^\w\-]/, "_")
          [sample, it]
        }
      ch_fasta_for_pangolin = CONSENSUS.out | mix(ch_extra_fasta_with_sample_name)
    } else {
      ch_fasta_for_pangolin = CONSENSUS.out
    }
    ch_fasta_for_pangolin \
      | PREPARE_FASTA_FOR_PANGOLIN \
      | collectFile(name: "sequences.fasta", newLine: true) \
      | PANGOLIN \
      | PANGOLIN_SUMMARY_FOR_MULTIQC
    ch_pangolin = PANGOLIN.out
    ch_pangolin_mqc = PANGOLIN_SUMMARY_FOR_MULTIQC.out
  } else {
    ch_pangolin = Channel.from([])
    ch_pangolin_mqc = Channel.from([])
  }

  if (params.tree) {
    include { MAFFT_MSA } from './processes/msa'
    include { IQTREE } from './processes/iqtree'
    include { BASIC_TREE_PLOT } from './processes/plots'
    // TODO: only run tree when number of samples is at least 2; IQ-TREE will only run if at least 3 sequences are provided as input
    ch_consensi = CONSENSUS.out | map { it[1] }
    if (params.tree_extra_fasta) {
      ch_extra_fasta = Channel.fromPath(params.tree_extra_fasta).splitFasta(file: true)
      ch_fasta_for_mafft = ch_consensi | mix(ch_extra_fasta) | collect
    } else {
      ch_fasta_for_mafft = ch_consensi | collect
    }
     
    // TODO: merge other sequences provided by user
    MAFFT_MSA(ch_fasta_for_mafft, ch_fasta)
    IQTREE(MAFFT_MSA.out, Channel.value(params.iqtree_model))
    // TODO: shiptv
    // TODO: merge optionally specified
    BASIC_TREE_PLOT(IQTREE.out, ch_pangolin)
    ch_basic_tree_plot = BASIC_TREE_PLOT.out
  } else {
    ch_basic_tree_plot = Channel.from([])
  }

  // Metagenomic classification by Kraken2
  if (params.kraken2_db && !params.skip_kraken) {
    KRAKEN2(file(params.kraken2_db), ch_reads)
  }

  Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt style=\"width:240px !important;\">$k</dt><dd style=\"margin-left:260px !important;\"><samp>${v != null ? v : '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'peterk87-nf-virontus-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'peterk87/nf-virontus Workflow Summary'
    section_href: 'https://github.com/peterk87/nf-virontus'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }
  SOFTWARE_VERSIONS()
  ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yaml")
  MULTIQC(
   ch_multiqc_config,
   MAP.out.stats.collect().ifEmpty([]),
   MOSDEPTH_GENOME.out.mqc.collect().ifEmpty([]),
   BCFTOOLS_STATS_POST_FILTER.out.collect().ifEmpty([]),
   ch_snpeff.ifEmpty([]),
   CONSENSUS_TO_MULTIQC_HTML.out.collect().ifEmpty([]),
   ch_pangolin_mqc.ifEmpty([]),
   ch_basic_tree_plot.ifEmpty([]),
   SOFTWARE_VERSIONS.out.software_versions_yaml.collect(),
   ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")
  )
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
