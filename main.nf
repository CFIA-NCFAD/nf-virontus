#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

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

  ${c_bul}De Novo Assembly Options:${c_reset}
    --do_unicycler_assembly       Assemble filtered reads using Unicycler? (default: ${params.do_unicycler_assembly})

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
summary['Ref Sequences FASTA'] = params.ref_fasta
if(params.bedfile) {
  summary['Primer Scheme'] = params.bedfile
}
summary['Consensus No Coverage'] = "<=${params.no_coverage}X positions replaced with '${params.no_cov_char}'"
summary['Consensus Low Coverage'] = "<${params.low_coverage}X positions replaced with '${params.low_cov_char}'"
summary['Centrifuge DB'] = params.centrifuge_db
summary['Kraken2 DB']   = params.kraken2_db
summary['Taxids'] = "Filtering for taxids belonging to $taxids"
summary['Unicycler Assembly?'] = params.do_unicycler_assembly ? "Yes" : "No"
if(params.do_unicycler_assembly) {
  summary['Unicycler Mode'] = params.unicycler_mode
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
log.info summary.collect { k,v -> "${k.padRight(22)}: $v" }.join("\n")
log.info "========================================="

//=============================================================================
// PROCESSES
//=============================================================================

// FASTA record to file
process REC2FASTA {
  tag "${record.id} - ${record.desc} - ${record.sequence.size()}"
  publishDir "${params.outdir}/refs", pattern: "*.fa", mode: 'copy'

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

process MAP {
  tag "$sample VS $ref_name"
  publishDir "${params.outdir}/mapping/$sample/bamfiles", pattern: "*.bam"

  input:
    tuple sample, 
          path(fastq), 
          path(ref_fasta)
  output:
    tuple sample, 
          path(ref_fasta),
          path(bam)

  script:
  ref_name = ref_fasta.getBaseName()
  bam = "${sample}-${ref_name}.bam"
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

process IVAR_TRIM {
  publishDir "${params.outdir}/mapping/$sample/bamfiles", pattern: "*.trim.bam"
  input:
    path(bedfile)
    tuple sample,
          path(ref_fasta),
          path(bam)
  output:
    tuple sample,
          path(ref_fasta),
          path(trimmed_bam)

  script:
  ref_name = ref_fasta.getBaseName()
  trimmed_bam = "${sample}-${ref_name}.trim.bam"
  """
  ivar trim \\
    -i $bam \\
    -b $bedfile \\
    -p trim -q 1 -m 20 -s 4 -e
  samtools sort -o $trimmed_bam trim.bam
  rm trim.bam
  """
}

process MAP_STATS {
  tag "$sample VS $ref_name"
  publishDir "${params.outdir}/mapping/$sample", mode: 'copy', pattern: "*.{tsv,flagstat,idxstats}"

  input:
    tuple val(sample),
          path(ref_fasta),
          path(bam)

  output:
    tuple val(sample),
          path(ref_fasta),
          path(bam),
          path(depths),
          path(flagstat),
          path(idxstats)
  script:
  ref_name = ref_fasta.getBaseName()
  depths = "${sample}-${ref_name}-depths.tsv"
  flagstat = "${sample}-${ref_name}.flagstat"
  idxstats = "${sample}-${ref_name}.idxstats"
  """
  samtools flagstat $bam > $flagstat
  samtools depth -a -d 0 $bam | perl -ne 'chomp \$_; print "${sample}\t\$_\n"' > $depths
  samtools idxstats $bam | head -n1 | perl -ne 'chomp \$_; print "${sample}\t\$_\n"' > $idxstats
  """
}

process MEDAKA {
  tag "$sample - $ref_name"
  publishDir "${params.outdir}/vcf", mode: 'copy', pattern: '*.vcf'

  input:
    tuple val(sample),
          path(ref_fasta),
          path(bam),
          path(depths),
          path(flagstat),
          path(idxstats)
  output:
    tuple val(sample),
          path(ref_fasta),
          path(bam),
          path(depths),
          path(vcf)
  script:
  ref_name = ref_fasta.getBaseName()
  vcf = "${sample}-${ref_name}.medaka.vcf"
  """
  samtools index $bam
  medaka consensus --chunk_len 800 --chunk_ovlp 400 $bam ${bam}.hdf
  medaka variant $ref_fasta ${bam}.hdf $vcf
  """
}

process LONGSHOT {
  tag "$sample - $ref_name"
  publishDir "${params.outdir}/vcf", mode: 'copy', pattern: '*.vcf'

  input:
    tuple val(sample),
          path(ref_fasta),
          path(bam),
          path(depths),
          path(medaka_vcf)
  output:
    tuple val(sample),
          path(ref_fasta),
          path(bam),
          path(depths),
          path(longshot_vcf)
  script:
  ref_name = ref_fasta.getBaseName()
  longshot_vcf = "${sample}-${ref_name}.longshot.vcf"
  script:
  """
  samtools faidx $ref_fasta
  samtools index $bam
  longshot -P 0 -F -A --no_haps \\
    --potential_variants $medaka_vcf \\
    --bam $bam \\
    --ref $ref_fasta \\
    --out $longshot_vcf
  """
}


// Filter for ALT allele variants that have greater depth than REF and
// that have greater than 2X coverage 
process BCF_FILTER {
  tag "$sample - $ref_name"
  publishDir "${params.outdir}/vcf",
    pattern: "*.filt.vcf",
    mode: 'copy'

  input:
    tuple val(sample),
          path(ref_fasta),
          path(bam),
          path(depths),
          path(vcf)
  output:
    tuple val(sample),
          path(ref_fasta),
          path(bam),
          path(depths),
          path(filt_vcf)
  script:
  ref_name = ref_fasta.getBaseName()
  filt_vcf = "${file(vcf).getBaseName()}.filt.vcf"
  """
  bcftools filter \\
    -e 'AC[0] >= AC[1] || AC[1]<=2' \\
    $vcf \\
    -Ov \\
    -o $filt_vcf
  """
}

process CONSENSUS {
  tag "$sample - $ref_name"
  publishDir "${params.outdir}/consensus", 
    pattern: "*.consensus.fasta",
    mode: 'copy'

  input:
    tuple val(sample),
          path(ref_fasta),
          path(bam),
          path(depths),
          path(filt_vcf)
  output:
    tuple val(sample),
          path(ref_fasta),
          path(bam),
          path(depths),
          path(consensus)

  script:
  ref_name = ref_fasta.getBaseName()
  consensus = "${sample}-${ref_name}.consensus.fasta"
  """
  vcf_consensus_builder \\
    -v $filt_vcf \\
    -d $depths \\
    -r $ref_fasta \\
    -o $consensus \\
    --low-coverage $params.low_coverage \\
    --no-coverage $params.no_coverage \\
    --low-cov-char $params.low_cov_char \\
    --no-cov-char $params.no_cov_char \\
    --sample-name $sample
  """
}

process KRAKEN2 {
  tag "$sample"
  publishDir "${params.outdir}/kraken2/results",
    pattern: "*-kraken2_results.tsv",
    mode: 'copy'
  publishDir "${params.outdir}/kraken2/reports",
    pattern: "*-kraken2_report.tsv",
    mode: 'copy'

  input:
    path(db)
    tuple val(sample), 
          path(reads)
  output:
    tuple val(sample),
          path(reads),
          path(results),
          path(report)

  script:
  results = "${sample}-kraken2_results.tsv"
  report = "${sample}-kraken2_report.tsv"
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
  tag "$sample"
  publishDir "${params.outdir}/centrifuge/$sample",
    pattern: "*.tsv",
    mode: 'copy'
  publishDir "${params.outdir}/centrifuge",
    pattern: "*-kreport.tsv",
    mode: 'copy'

  input:
    tuple db_name, 
          path(db)
    tuple val(sample),
          path(reads)
  output:
    tuple val(sample),
          path(reads),
          path(results),
          path(kreport)

  script:
  results = "${sample}-centrifuge_results.tsv"
  kreport = "${sample}-kreport.tsv"
  """
  centrifuge \\
    -x ${db}/${db_name} \\
    -U $reads \\
    -S $results \\
    -p ${task.cpus}
  centrifuge-kreport -x ${db}/${db_name} $results > $kreport
  """
}

process FILTER_READS_BY_CLASSIFICATIONS {
  tag "$sample"
  publishDir "${params.outdir}/filtered_reads/", 
    pattern: "*.viral_unclassified.fastq", 
    mode: 'copy'
  errorStrategy 'ignore'

  input:
    tuple sample,
          path(reads),
          path(kraken2_results),
          path(kraken2_report),
          path(centrifuge_results),
          path(centrifuge_report)
  output:
    tuple sample,
          path(filtered_reads) optional true

  script:
  filtered_reads = "${sample}.viral_unclassified.fastq"
  taxids_arg = taxids ? " --taxids $taxids" : ""
  """
  filter_classified_reads \\
    ${taxids_arg} \\
    -i $reads \\
    -o $filtered_reads \\
    -c $centrifuge_results \\
    -C $centrifuge_report \\
    -k $kraken2_results \\
    -K $kraken2_report
  """
}

process UNICYCLER_ASSEMBLY {
  tag "$sample"
  publishDir "${params.outdir}/assemblies/unicycler/$sample", mode: 'copy'
  errorStrategy 'ignore'

  input:
    tuple val(sample),
          path(reads)
  output:
    tuple val(sample),
          val('unicycler'),
          path("${sample}/")

  script:
  output_contigs = "${sample}-assembly.fasta"
  output_gfa = "${sample}-assembly.gfa"
  output_unicycler_log = "${sample}-unicycler.log"
  """
  unicycler -t ${task.cpus} --mode ${params.unicycler_mode} -o $sample -l $reads
  ln -s ${sample}/assembly.fasta $output_contigs
  ln -s ${sample}/assembly.gfa $output_gfa
  ln -s ${sample}/unicycler.log $output_unicycler_log
  """
}

//=============================================================================
// WORKFLOW
//=============================================================================
workflow {

  // Reference genome FASTA input channel
  Channel.fromPath( params.ref_fasta )
    .splitFasta( record: [id: true, desc: true, sequence: true] ) \
    | REC2FASTA 

  // Map each sample's reads against each reference genome sequence
  // - Combine each ref seq with each sample's reads
  // - map reads against ref
  Channel.fromPath(params.reads)
    .map { [file(it).getBaseName(), it] }
    .set { ch_reads }

  ch_reads | combine(REC2FASTA.out) | MAP 

  // Trim primer sequences from read alignments if primer scheme BED file provided 
  if (params.bedfile) {
    Channel.value( file(params.bedfile) )
      .set { ch_bedfile}
    IVAR_TRIM(ch_bedfile, MAP.out) | MAP_STATS
  } else {
    MAP_STATS(MAP.out)
  }

  MAP_STATS.out \
    | filter { 
      // Filter for alignments that did have some reads mapping to the ref genome
      depth_linecount = file(it[3]).readLines().size()
      if (depth_linecount == 1) {
        println "No reads from \"${it[0]}\" mapped to reference ${it[1]}"
      }
      depth_linecount > 2
    } \
    | MEDAKA \
    | LONGSHOT \
    | BCF_FILTER \
    | CONSENSUS



  if (params.kraken2_db) {
    // Kraken2 DB input channel
    Channel.value( file(params.kraken2_db) )
      .set { ch_kraken2_db }

    // Metagenomic classification by Kraken2 and Centrifuge
    KRAKEN2(ch_kraken2_db, ch_reads)  
  }
  
  if (params.centrifuge_db) {
    // Centrifuge DB input channel
    Channel.value( 
      [ 
        file(params.centrifuge_db).getName(), 
        file(params.centrifuge_db).getParent() 
      ] )
      .set { ch_centrifuge_db }

    CENTRIFUGE(ch_centrifuge_db, ch_reads)
  }
  
  if (params.kraken2_db && params.centrifuge_db) {
    // Join Kraken2 and Centrifuge classification results
    ch_k2_cent_res = KRAKEN2.out.join(CENTRIFUGE.out, remainder: true)
      .map { sample, reads, kraken2_results, kraken2_report, _r, centrifuge_results, centrifuge_kreport -> 
        [
          sample, 
          reads, 
          kraken2_results, 
          kraken2_report, 
          centrifuge_results, 
          centrifuge_kreport
        ]
      }

    FILTER_READS_BY_CLASSIFICATIONS(ch_k2_cent_res)
    if (params.do_unicycler_assembly) {
      FILTER_READS_BY_CLASSIFICATIONS.out | UNICYCLER_ASSEMBLY
    }
  } else {
    if (params.do_unicycler_assembly) {
      ch_reads | UNICYCLER_ASSEMBLY
    }
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
