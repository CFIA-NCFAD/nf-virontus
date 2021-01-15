#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
  helpMessage;
  checkFileExists;
  checkKraken2Db;
  checkTaxids;
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
  MAP
} from './processes/mapping'
include {
  IVAR_TRIM
} from './processes/ivar'
include {
  MEDAKA_LONGSHOT;
  VARIANT_FILTER
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

//=============================================================================
// WORKFLOW
//=============================================================================
workflow {
  // taxids = checkTaxids(params.taxids)
  // if (params.centrifuge_db) checkCentrifugeDb(params.centrifuge_db)
  // if (params.kraken2_db) checkKraken2Db(params.kraken2_db)

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

  // Map each sample's reads against each reference genome sequence
  // - Combine each ref seq with each sample's reads
  // - map reads against ref
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
  ch_fasta = Channel.value(file(fasta))
  if (gff) {
    MAKE_SNPEFF_DB(Channel.value([index_base, file(fasta), file(gff)]))
  }

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

  ch_bam | MEDAKA_LONGSHOT
  SNPEFF(MEDAKA_LONGSHOT.out, MAKE_SNPEFF_DB.out)
  VARIANT_FILTER(MEDAKA_LONGSHOT.out, params.major_allele_fraction)

  VARIANT_FILTER.out | join(ch_depths) | (CONSENSUS & COVERAGE_PLOT)

  // // Metagenomic classification by Kraken2
  // if (params.kraken2_db) {
  //   KRAKEN2(file(params.kraken2_db), ch_reads)  
  // }
  
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
