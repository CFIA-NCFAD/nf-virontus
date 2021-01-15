process KRAKEN2 {
  tag "$sample"
  publishDir "${params.outdir}/kraken2/results",
    pattern: "*-kraken2_results.tsv",
    mode: params.publish_dir_mode
  publishDir "${params.outdir}/kraken2/reports",
    pattern: "*-kraken2_report.tsv",
    mode: params.publish_dir_mode

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