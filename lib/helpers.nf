/*
 * Helper functions for peterk87/nf-virontus
 */

// TODO: update help
def helpMessage() {
  // Log colors ANSI codes
  c_reset = params.monochrome_logs ? '' : "\033[0m";
  c_bold = params.monochrome_logs ? '' : "\033[1m";
  c_dim = params.monochrome_logs ? '' : "\033[2m";
  c_block = params.monochrome_logs ? '' : "\033[3m";
  c_ul = params.monochrome_logs ? '' : "\033[4m";
  c_black = params.monochrome_logs ? '' : "\033[0;30m";
  c_red = params.monochrome_logs ? '' : "\033[0;31m";
  c_green = params.monochrome_logs ? '' : "\033[0;32m";
  c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
  c_blue = params.monochrome_logs ? '' : "\033[0;34m";
  c_purple = params.monochrome_logs ? '' : "\033[0;35m";
  c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
  c_white = params.monochrome_logs ? '' : "\033[0;37m";
  c_bul = c_bold + c_ul;
  is_viruses = (params.taxids == 10239) ? " (Viruses)" : ""
  log.info"""
  =${c_dim}=================================================================${c_reset}
  ${c_blue+c_bold}${workflow.manifest.name}${c_reset}   ~  version ${c_purple}${workflow.manifest.version}${c_reset}
  ${c_dim}==================================================================${c_reset}

    ${c_ul}Git info:${c_reset} $workflow.repository - $workflow.revision [$workflow.commitId]

  ${c_bul}Usage:${c_reset}
  Given some barcoded and demultiplexed reads, the typical command for running the pipeline is as follows:
  
    nextflow run ${workflow.manifest.name} \\
      ${c_red}--reads "${params.reads}"${c_reset} \\
      ${c_green}--outdir ${params.outdir}${c_reset} \\
      --ref_fasta refs.fa \\
      -profile singularity # recommended to run with Singularity

  The above ${c_bul}assumes${c_reset} that you have a ${c_cyan}Centrifuge DB${c_reset} and ${c_purple}Kraken2 DB${c_reset} located at
  ${c_cyan}/opt/DB/centrifuge/nt-2018-03-03/nt${c_reset} and ${c_purple}/opt/DB/kraken2/standard2${c_reset}, 
  respectively, ${c_bul}OR${c_reset} that you have set ${c_cyan}\$CENTRIFUGE_DB${c_reset} and ${c_purple}\$KRAKEN2_DB${c_reset} env 
  variables. It also assumes that you have ${c_yellow+c_bul}Singularity${c_reset} installed on your
  local machine and will automatically pull and use the Singularity image for
  this workflow from Singularity-Hub.org.

  ${c_yellow+c_bold+c_block}NOTE:${c_yellow} For best results, please ensure you have ${c_bul}Singularity${c_yellow} installed prior to running this workflow.${c_dim}(https://sylabs.io/guides/3.3/user-guide/quick_start.html#quick-installation-steps)${c_reset}

  Note: 
  The argument supplied to "--reads" must be quoted if using "*" and other 
  characters and symbols that could be shell expanded!

  ${c_bul}Mandatory Options:${c_reset}
    ${c_red}--reads${c_reset}   Input reads directory and pattern (default: ${c_red}"${params.reads}"${c_reset})
    --ref_fasta      Reference genomes multiFASTA file (one or more references
                     in a single file) (default: "${file(params.ref_fasta)}")
  ${c_bul}Amplicon Sequencing Options:${c_reset}
    --bedfile        BED format file with amplicon sequencing primers info (optional). 
                     Produced as output from PrimalScheme.
  ${c_bul}Consensus Generation Options:${c_reset}
    --low_coverage   Low coverage threshold (default=${params.low_coverage}).
                     Replace consensus sequence positions below this depth
                     threshold with a low coverage character 
                     (see ${c_dim}--low_cov_char${c_reset})
    --no_coverage    No coverage threshold (default=${params.no_coverage}).
                     Replace consensus sequence positions with less than or 
                     equal this depth with a no coverage character 
                     (see ${c_dim}--no_cov_char${c_reset})
    --low_cov_char   Low coverage character (default="${params.low_cov_char}")
    --no_cov_char    No coverage character (default="${params.no_cov_char}")

  SARS-CoV-2:
  ${params.genomes[params.genome]}

  ${c_bul}Taxonomic Classification Options:${c_reset}
    ${c_cyan}--centrifuge_db${c_reset}   Path to Centrifuge DB and prefix. If not specified, will 
                      try to get from \$CENTRIFUGE_DB env variable or see if
                      "/opt/DB/centrifuge/nt-2018-03-03/nt" exists.
                      (default: ${c_cyan}${params.centrifuge_db}${c_reset})
    ${c_purple}--kraken2_db${c_reset}      Path to Kraken2 DB directory. . If not specified, will 
                      try to get from \$KRAKEN2_DB env variable or see if
                      "/opt/DB/kraken2/standard2" exists.
                      (default: ${c_purple}${params.kraken2_db}${c_reset})
    --taxids          Taxonomic IDs to filter reads by. Multiple taxids should
                      be delimited by commas (`--taxids 1,2,3`). To disable 
                      filtering of reads based on taxids, do not provide a
                      value for the `--taxids` argument:
                      `nextflow run ... --taxids --reads ...`
                      (default: ${params.taxids}${is_viruses})
    --exclude_unclassified_reads  Exclude unclassified reads from taxonomic
                                  classification filtered reads (default: false)

  ${c_bul}De Novo Assembly Options:${c_reset}
    --do_unicycler_assembly       Assemble filtered reads using Unicycler? (default: ${params.do_unicycler_assembly})

  ${c_bul}Cluster Options:${c_reset}
    --slurm_queue     Name of SLURM queue to run workflow on; use with ${c_dim}-profile slurm${c_reset}

  ${c_bul}Other Options:${c_reset}
    ${c_green}--outdir${c_reset}          The output directory where the results will be saved
                      (default: ${c_green}${params.outdir}${c_reset})
    -w/--work-dir     The temporary directory where intermediate data will be 
                      saved (default: ${workflow.workDir})
    -profile          Configuration profile to use. [standard, singularity, 
                      conda, slurm] (default '${workflow.profile}')
    --tracedir        Pipeline run info output directory (default: 
                      ${params.tracedir})

  Note: 
  It is recommended that this workflow be executed with Singularity using the 
  Singularity profile (`-profile singularity`) for maximum reproducibility and
  ease of execution on different platforms.
  """.stripIndent()
}

//=============================================================================
// User input validation helper functions
//=============================================================================

def checkFileExists(file_path) {
  f = file(file_path)
  if ( !f.isFile() || !f.exists() ) {
    exit 1, "File '$file_path' does not exist!"
  }
}

def checkCentrifugeDb(centrifuge_db) {
  file_centrifuge_db = file(centrifuge_db)
  prefix = file_centrifuge_db.getName()
  centrifuge_dir = file_centrifuge_db.getParent()
  if ( !centrifuge_dir.isDirectory() || !centrifuge_dir.exists() ) {
    exit 1, "Centrifuge DB does not exist at '$centrifuge_dir'! Please specify a valid Centrifuge DB."
  }
  any_valid = false
  centrifuge_dir.eachFile { f ->
    if ( f.isFile() ) {
      if ( f.getName() =~ /^$prefix/ && f.getExtension() == 'cf') {
        any_valid = true
      }
    }
  }
  if ( !any_valid ) {
    exit 1, "No valid Centrifuge DB files with prefix '$prefix' in '$centrifuge_dir' and extension 'cf'! Please specify a valid Centrifuge classification DB directory and prefix."
  }
}

def checkKraken2Db(kraken2_db) {
  kraken2_db_dir = file(kraken2_db)
  if ( !kraken2_db_dir.isDirectory() ) {
    exit 1, "The Kraken2 DB must be a directory! '$kraken2_db' is not a directory!"
  }
  if ( !kraken2_db_dir.exists() ) {
    exit 1, "The Kraken2 DB must be an existing directory! '$kraken2_db' does not exist!"
  }
}

// Check that all taxids are integers delimited by commas
def checkTaxids(taxids) {
  if (taxids instanceof Boolean || taxids.toString().isEmpty()) {
    return null
  } else if (taxids.toString().isInteger()) {
    return taxids.toString()
  } else {
    taxids_list = taxids.toString()
      .split(',')
      .collect { it.strip() }
      .findAll { it != '' }
    if (!taxids_list.every { it.isInteger() }) {
      exit 1, "Not every element in `--taxids` is an integer!"
    }
    return taxids_list.join(',')
  }
}

def check_sample_sheet(LinkedHashMap sample_sheet) {
  // Check that each entry from a sample sheet
  reads_file = file(sample_sheet.reads, checkIfExists: true)
  reads = sample_sheet.reads ? reads_file : null
  if (reads == null) {
    exit 1, "The Nanopore reads FASTQ file or directory specified for ${sample_sheet.sample} does not exist! Please check the sample sheet '$params.sample_sheet'"
  }
  if (reads_file.isDirectory()) {
    fqs = []
    reads_file.eachFile {
      fname = it.getName()
      if (fname ==~ /.*\.fastq(\.gz)?$/) {
        fqs << it
      }
    }
    if (fqs.size() == 0) {
      exit 1, "Sample '${sample_sheet.sample}' input specified as directory '$reads_file' with no FASTQ files! Please check the sample sheet '$params.sample_sheet'\n${fqs}\n${reads_file.listFiles()}"
    }
    reads = fqs
  }
  return [ sample_sheet.sample, reads ]
}
