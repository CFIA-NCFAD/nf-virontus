process KRAKEN2 {
  tag "$sample"
  label "process_high"
  publishDir "${params.outdir}/kraken2/",
    mode: params.publish_dir_mode

  input:
  path(db)
  tuple val(sample), path(reads)
  
  output:
  tuple val(sample), path('*.unclassified.fastq.gz'), emit: unclassified_reads
  tuple val(sample), path('*.classified.fastq.gz'), emit: classified_reads
  tuple val(sample), path(report), emit: report
  tuple val(sample), path(results), emit: results

  script:
  results = "${sample}.kraken2.results.tsv"
  report = "${sample}.kraken2.report.txt"
  unclassified_reads = "${sample}.unclassified.fastq"
  classified_reads = "${sample}.classified.fastq"
  """
  kraken2 \\
    --db $db \\
    --unclassified-out $unclassified_reads \\
    --classified-out $classified_reads \\
    --report $report \\
    --output $results \\
    --threads ${task.cpus} \\
    --gzip-compressed \\
    $reads
  
  pigz *.fastq
  """
}
