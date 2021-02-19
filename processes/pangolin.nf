process PREPARE_FASTA_FOR_PANGOLIN {
  input:
  tuple val(sample), path(fasta)

  output:
  path(output)

  script:
  output = "${sample}.fasta"
  """
  cat $fasta > $output
  sed -i 's/^>.*/>${sample}/g' $output
  """
}

process PANGOLIN {
  container "covlineages/pangolin:v2.3.0"

  publishDir "${params.outdir}/pangolin",
             pattern: "*.csv",
             mode: params.publish_dir_mode

  publishDir "${params.outdir}/pangolin/logs",
             pattern: ".command.log",
             saveAs: { "pangolin.log" },
             mode: params.publish_dir_mode

  input:
  path(fasta)

  output:
  path(output)

  script:
  output = "pangolin.lineage_report.csv"
  """
  # update pangolin lineages in case of updates
  git clone --depth 1 https://github.com/cov-lineages/pangoLEARN.git
  pip install --upgrade pangoLEARN/

  pangolin \\
    --verbose \\
    $fasta \\
    --outfile $output
  """
}

process PANGOLIN_SUMMARY_FOR_MULTIQC {
  input:
  path(input)

  output:
  path("pangolin_mqc.txt")

  script:
  """
  summarize_pangolin_results.py ./ pangolin_mqc.txt
  """
}
