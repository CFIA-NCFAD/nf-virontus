process CONSENSUS {
  tag "$sample - $ref_name"
  publishDir "${params.outdir}/consensus", 
    pattern: "*.consensus.fasta",
    mode: params.publish_dir_mode

  input:
  tuple val(sample),
        path(ref_fasta),
        path(vcf),
        path(depths)
  output:
  path(consensus)

  script:
  ref_name = ref_fasta.getBaseName()
  consensus = "${sample}-${ref_name}.consensus.fasta"
  """
  vcf_consensus_builder \\
    -v $vcf \\
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
