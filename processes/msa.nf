process MAFFT_MSA {
  publishDir "${params.outdir}/msa",
             pattern: "mafft.fasta",
             mode: params.publish_dir_mode
  input:
  path(fastas)
  path(ref_fasta)

  output:
  path('mafft.fasta')

  script:
  """
  # ensure that ref fasta is first
  cat $ref_fasta > input.fasta
  # only keep seq id
  cat $fastas | sed -r 's/^>(\\S+)\\s.*/>\\1/g' >> input.fasta
  mafft \\
    --thread ${task.cpus} \\
    --auto \\
    input.fasta > mafft.fasta
  """
}
