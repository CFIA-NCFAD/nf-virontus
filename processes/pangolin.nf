process PANGOLIN {
  tag "$sample"
  // TODO: keep up-to-date with new releases of Pangolin
  container "covlineages/pangolin:v2.1.7"

  publishDir "${params.outdir}/pangolin",
             pattern: "*.csv",
             mode: params.publish_dir_mode
  publishDir "${params.outdir}/pangolin/logs",
             pattern: ".command.log",
             saveAs: { "pangolin-${sample}.log" },
             mode: params.publish_dir_mode

  input:
  tuple val(sample), path(fasta)

  output:
  tuple val(sample), path(output)

  script:
  output = "${sample}.pangolin.lineage_report.csv"
  """
  cp $fasta tmp.fasta
  sed -i 's/^>.*/>${sample}/g' tmp.fasta
  pangolin \\
    --verbose \\
    -t ${task.cpus} \\
    tmp.fasta \\
    --outfile $output
  """
}

process PANGOLIN_SUMMARY_FOR_MULTIQC {
  input:
  tuple val(value), path(input)

  output:
  path("pangolin_mqc.txt")

  script:
  """
  summarize_pangolin_results.py ./ pangolin_mqc.txt
  """
}
