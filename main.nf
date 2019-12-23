#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

min_nextflow_version = '19.10'
if (!nextflow.version.matches("${min_nextflow_version}+")) {
  log.error "This workflow requires Nextflow version $min_nextflow_version or greater. You are running version $nextflow.version! Please install a newer version or run `nextflow self-update` to self-update Nextflow."
  exit 1
}

valid_barcode_kits = ['EXP-NBD103', 'EXP-NBD104', 'EXP-NBD114', 'EXP-PBC001',
  'EXP-PBC096', 'SQK-16S024', 'SQK-LWB001', 'SQK-PBK004', 'SQK-RAB201',
  'SQK-RAB204', 'SQK-RBK001', 'SQK-RBK004', 'SQK-RLB001', 'SQK-RPB004',
  'VSK-VMK001', 'VSK-VMK002'] as Set

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
  The typical command for running the pipeline is as follows:
  
    nextflow run ${workflow.manifest.name} \\
      ${c_red}--reads "${params.reads}"${c_reset} \\
      ${c_green}--outdir ${params.outdir}${c_reset} \\
      --ref_fasta refs.fa \\
      --barcode_kit $params.barcode_kit \\
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

  ${c_bul}Barcoding Options:${c_reset}
    --barcode_kit    Nanopore barcoding kit (default: "$params.barcode_kit")

  ${c_bul}Cluster Options:${c_reset}
    --slurm_queue     Name of SLURM queue to run workflow on; use with ${c_dim}-profile slurm${c_reset}

    
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

  Nanopore barcoding kit (`--barcode_kit`) must be one of:
  ${valid_barcode_kits.join(', ')}
  """.stripIndent()
}
//=============================================================================
// Help info
//=============================================================================
// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
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

//=============================================================================
// Check user input params
//=============================================================================
if (workflow.profile == 'slurm' && params.slurm_queue == "") {
  log.error "You must specify a valid SLURM queue (e.g. '--slurm_queue <queue name>' (see `\$ sinfo` output for available queues)) to run this workflow with the 'slurm' profile!"
  exit 1
}

if (!(params.barcode_kit in valid_barcode_kits)) {
  log.error "Invalid barcode kit '${params.barcode_kit}'! Must be one of '${valid_barcode_kits}'"
  exit 1
}

checkFileExists(params.ref_fasta)
taxids = checkTaxids(params.taxids)
if (params.centrifuge_db) checkCentrifugeDb(params.centrifuge_db)
if (params.kraken2_db) checkKraken2Db(params.kraken2_db)

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

//=============================================================================
// LOG EXECUTION START PARAMS
//=============================================================================
log.info """=======================================================
${workflow.manifest.name} v${workflow.manifest.version}
======================================================="""
def summary = [:]
summary['Pipeline Name']  = workflow.manifest.name
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Reads']        = params.reads
summary['Reference Genomes FASTA'] = params.ref_fasta
summary['Centrifuge DB'] = params.centrifuge_db
summary['Kraken2 DB']   = params.kraken2_db
summary['Taxids'] = taxids
summary['Assembly with Unicycler'] = params.do_unicycler_assembly
if(params.do_unicycler_assembly) {
  summary['Unicycler Mode'] = params.unicycler_mode
}
if(params.do_flye_assembly) {
  summary['Assembly with Flye'] = params.do_flye_assembly
}
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
summary['Command-Line']   = workflow.commandLine
summary['Nextflow version'] = workflow.nextflow.version
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

//=============================================================================
// PROCESSES
//=============================================================================

// FASTA record to file
process REC2FASTA {
  tag "${record.id} - ${record.desc} - ${record.sequence.size()}"

  input:
    val(record)
  output:
    path(fasta)

  script:
  fasta = "${record.id}.fa"
  """
  cat > $fasta << EOF
>${record.id} ${record.desc}
${record.sequence}
EOF
  """
}

process BARCODER {
  tag "$fastq_name"

  input:
    file(fastq)

  output:
    path('**/*.fastq.gz')

  script:
  fastq_name = fastq.getBaseName()
  """
  guppy_barcoder \\
    -i ./ \\
    -s ./ \\
    --barcode_kits ${params.barcode_kit} \\
    --trim_barcodes \\
    --detect_mid_strand_barcodes \\
    --compress_fastq \\
    -q 0 \\
    -t ${task.cpus}
  """
}

process CAT_FASTQS {
  tag "$barcode"

  input:
    tuple val(barcode),
          path('reads*.fastq.gz')
  output:
    tuple val(barcode),
          path("${barcode}.fastq.gz")

  """
  cat reads*.fastq.gz > ${barcode}.fastq.gz
  """
}

process MAP {
  tag "$barcode VS $ref_fasta"
  publishDir "${params.outdir}/mapping/$barcode/bamfiles", pattern: "*.bam"
  publishDir "${params.outdir}/refs", pattern: "*.fasta", mode: 'copy'

  input:
    tuple barcode, 
          path(fastq), 
          path(ref_fasta)
  output:
    tuple barcode, 
          path(ref_fasta),
          path(bam)

  script:
  ref_name = ref_fasta.getBaseName()
  bam = "${barcode}-${ref_name}.bam"
  """
  minimap2 \\
    -ax map-ont \\
    -t${task.cpus} \\
    $ref_fasta \\
    $fastq \\
    | samtools sort -@${task.cpus} \\
    | samtools view -F4 -b -o $bam -
  """
}

process MAP_STATS {
  tag "$barcode VS segment $ref_name"
  publishDir "${params.outdir}/mapping/$barcode", mode: 'copy', pattern: "*.{tsv,flagstat,idxstats}"

  input:
    tuple val(barcode),
          path(ref_fasta),
          path(bam)

  output:
    tuple val(barcode),
          path(ref_fasta),
          path(bam),
          path(depths),
          path(flagstat),
          path(idxstats)
  script:
  ref_name = ref_fasta.getBaseName()
  depths = "${barcode}-${ref_name}-depths.tsv"
  flagstat = "${barcode}-${ref_name}.flagstat"
  idxstats = "${barcode}-${ref_name}.idxstats"
  """
  samtools flagstat $bam > $flagstat
  samtools depth -a -d 0 $bam | perl -ne 'chomp \$_; print "${barcode}\t\$_\n"' > $depths
  samtools idxstats $bam | head -n1 | perl -ne 'chomp \$_; print "${barcode}\t\$_\n"' > $idxstats
  """
}

//TODO: filter empty mpileup results or ignore errors
process BCF_CALL {
  tag "$barcode - $ref_name"
  publishDir "${params.outdir}/bcftools/call", mode: 'copy', pattern: '*.bcf'

  input:
    tuple val(barcode),
          path(ref_fasta),
          path(bam),
          path(depths),
          path(flagstat),
          path(idxstats)
  output:
    tuple val(barcode),
          path(ref_fasta),
          path(bam),
          path(depths),
          path(bcf)
  script:
  ref_name = ref_fasta.getBaseName()
  bcf = "${barcode}-${ref_name}.bcf"
  """
  bcftools mpileup --threads ${task.cpus} \\
      -f $ref_fasta \\
      $bam \\
      -Q 3 \\
      -Ou \\
  | bcftools call --threads ${task.cpus} \\
    -mv - \\
    -o $bcf
  """
}

// Filter for variants that meet the following criteria:
// - Maximum fraction of reads supporting an indel is 0.5 or greater (IMF)
// - Depth of coverage of ALT allele is greater than REF allele coverage
process BCF_FILTER {
  tag "$barcode - $ref_name"
  publishDir "${params.outdir}/bcftools/call",
    pattern: "*.filt.vcf.gz",
    mode: 'copy'

  input:
    tuple val(barcode),
          path(ref_fasta),
          path(bam),
          path(depths),
          path(bcf)
  output:
    tuple val(barcode),
          path(ref_fasta),
          path(bam),
          path(depths),
          path(filt_vcfgz)
  script:
  ref_name = ref_fasta.getBaseName()
  filt_vcfgz = "${barcode}-${ref_name}.filt.vcf.gz"
  """
  bcftools norm --threads ${task.cpus} \\
    -f $ref_fasta \\
    $bcf \\
    -Ob \\
  | bcftools filter \\
    -e 'IMF<0.5 || (DP4[0]+DP4[1])>(DP4[2]+DP4[3])' \\
    - \\
    -Oz \\
    -o $filt_vcfgz
  """
}

process CONSENSUS {
  tag "$barcode - $ref_name"
  publishDir "${params.outdir}/bcftools/consensus", 
    pattern: "*.consensus.fasta",
    mode: 'copy'

  input:
    tuple val(barcode),
          path(ref_fasta),
          path(bam),
          path(depths),
          path(filt_vcfgz)
  output:
    tuple val(barcode),
          path(ref_fasta),
          path(bam),
          path(depths),
          path(consensus)

  script:
  ref_name = ref_fasta.getBaseName()
  consensus = "${barcode}-${ref_name}.consensus.fasta"
  """
  bcftools index $filt_vcfgz
  bcftools consensus -f $ref_fasta $filt_vcfgz > $consensus
  """
}

process FIX_CONSENSUS {
  publishDir "${params.outdir}/consensus",
    pattern: "*.fixed_consensus.fasta",
    mode: 'copy'

  input:
    tuple val(barcode),
          path(ref_fasta),
          path(bam),
          path(depths),
          path(consensus)
  output:
    tuple val(barcode),
          path(ref_fasta),
          path(bam),
          path(depths),
          path(fixed_consensus) optional true

  script:
  ref_name = ref_fasta.getBaseName()
  fixed_consensus = "${barcode}-${ref_name}.fixed_consensus.fasta"
  low_cov_char = "N"
  cov_threshold = 1
  """
  fix_consensus.py \\
    -s $consensus \\
    -d $depths \\
    -o $fixed_consensus \\
    --cov-threshold $cov_threshold \\
    --low-cov-char $low_cov_char \\
    --sample-name "$barcode"
  """
}

process KRAKEN2 {
  tag "$barcode"
  publishDir "${params.outdir}/kraken2/results",
    pattern: "*-kraken2_results.tsv",
    mode: 'copy'
  publishDir "${params.outdir}/kraken2/reports",
    pattern: "*-kraken2_report.tsv",
    mode: 'copy'

  input:
    path(db)
    tuple val(barcode), 
          path(reads)
  output:
    tuple val(barcode),
          path(reads),
          path(results),
          path(report)

  script:
  results = "${barcode}-kraken2_results.tsv"
  report = "${barcode}-kraken2_report.tsv"
  """
  kraken2 \\
    --threads ${task.cpus} \\
    --memory-mapping \\
    --db ./${db}/ \\
    --report ${report} \\
    --output ${results} \\
    $reads
  """
}

process CENTRIFUGE {
  tag "$barcode"
  publishDir "${params.outdir}/centrifuge/$barcode",
    pattern: "*.tsv",
    mode: 'copy'
  publishDir "${params.outdir}/centrifuge",
    pattern: "*-kreport.tsv",
    mode: 'copy'

  input:
    tuple db_name, 
          path(db)
    tuple val(barcode),
          path(reads)
  output:
    tuple val(barcode),
          path(reads),
          path(results),
          path(kreport)

  script:
  results = "${barcode}-centrifuge_results.tsv"
  kreport = "${barcode}-kreport.tsv"
  """
  centrifuge \\
    -x ${db}/${db_name} \\
    -U $reads \\
    -S $results \\
    --mm
  centrifuge-kreport -x ${db}/${db_name} $results > $kreport
  """
}

process FILTER_READS_BY_CLASSIFICATIONS {
  tag "$barcode"
  publishDir "${params.outdir}/filtered_reads/", 
    pattern: "*.viral_unclassified.fastq", 
    mode: 'copy'

  input:
    tuple barcode,
          path(reads),
          path(kraken2_results),
          path(kraken2_report),
          path(centrifuge_results),
          path(centrifuge_report)
  output:
    tuple barcode,
          path(filtered_reads) optional true

  script:
  filtered_reads = "${barcode}.viral_unclassified.fastq"
  taxids_arg = taxids ? " --taxids $taxids" : ""
  """
  filter_classified_reads \\
    ${taxids_arg}
    -i $reads \\
    -o $filtered_reads \\
    -c $centrifuge_results \\
    -C $centrifuge_report \\
    -k $kraken2_results \\
    -K $kraken2_report
  """
}

process UNICYCLER_ASSEMBLY {
  tag "$barcode"
  publishDir "${params.outdir}/assemblies/unicycler/$barcode", mode: 'copy'
  errorStrategy 'ignore'

  input:
    tuple val(barcode),
          path(reads)
  output:
    tuple val(barcode),
          val('unicycler'),
          path("${barcode}/")

  script:
  output_contigs = "${barcode}-assembly.fasta"
  output_gfa = "${barcode}-assembly.gfa"
  output_unicycler_log = "${barcode}-unicycler.log"
  """
  unicycler -t ${task.cpus} --mode ${params.unicycler_mode} -o $barcode -l $reads
  ln -s ${barcode}/assembly.fasta $output_contigs
  ln -s ${barcode}/assembly.gfa $output_gfa
  ln -s ${barcode}/unicycler.log $output_unicycler_log
  """
}

// TODO: try flye assembly with different genome sizes derived from classification results?
process FLYE_ASSEMBLY {
  tag "$barcode"
  publishDir "${params.outdir}/assemblies/flye/$barcode", mode: 'copy'
  errorStrategy 'ignore'

  input:
    tuple val(barcode),
          path(reads)
  output:
    path("out/")

  script:
  """
  flye \\
    --nano-raw $reads \\
    --out-dir out \\
    --genome-size $params.expected_genome_size
  """
}

//=============================================================================
// WORKFLOW
//=============================================================================
workflow {

  // Centrifuge DB input channel
  Channel.value( 
    [ 
      file(params.centrifuge_db).getName(), 
      file(params.centrifuge_db).getParent() 
    ] )
    .set { ch_centrifuge_db }

  // Kraken2 DB input channel
  Channel.value( file(params.kraken2_db) )
    .set { ch_kraken2_db }

  // Reference genome FASTA input channel
  Channel.fromPath( params.ref_fasta )
    .splitFasta( record: [id: true, desc: true, sequence: true] ) \
    | REC2FASTA 

  Channel.fromPath(params.reads) \
    | BARCODER \
    | flatMap \
    | map { [file(it).getParent().getName(), it] } \
    | groupTuple \
    | CAT_FASTQS

  // Map barcoded reads against each reference genome sequence
  // - Combine each ref seq with each barcoded reads
  // - map reads against ref
  // - samtools stats for mapping
  // - BCF mpileup and calling
  // - BCF filtering
  // - bcftools consensus
  // - fix consensus (mask low coverage positions)
  CAT_FASTQS.out \
    | combine(REC2FASTA.out) \
    | MAP \
    | MAP_STATS \
    | filter { 
      // Filter for alignments that did have some reads mapping to the ref genome
      depth_linecount = file(it[3]).readLines().size()
      if (depth_linecount == 1) {
        println "No reads from \"${it[0]}\" mapped to reference ${it[1]}"
      }
      depth_linecount > 2
    } \
    | BCF_CALL \
    | BCF_FILTER \
    | CONSENSUS \
    | FIX_CONSENSUS

  // Metagenomic classification by Kraken2 and Centrifuge
  KRAKEN2(ch_kraken2_db, CAT_FASTQS.out)
  CENTRIFUGE(ch_centrifuge_db, CAT_FASTQS.out)
  // Join Kraken2 and Centrifuge classification results
  ch_k2_cent_res = KRAKEN2.out.join(CENTRIFUGE.out, remainder: true)
    .map { barcode, reads, kraken2_results, kraken2_report, _r, centrifuge_results, centrifuge_kreport -> 
      [
        barcode, 
        reads, 
        kraken2_results, 
        kraken2_report, 
        centrifuge_results, 
        centrifuge_kreport
      ]
    }

  FILTER_READS_BY_CLASSIFICATIONS(ch_k2_cent_res)
  if (params.do_unicycler_assembly) {
    UNICYCLER_ASSEMBLY(FILTER_READS_BY_CLASSIFICATIONS.out)
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
