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
  # keep only the seq id; remove all other text after fasta seq id
  cat $fastas | sed -r 's/^>(\\S+)\\s.*/>\\1/g' > input.fasta
  mafft \\
    --thread ${task.cpus} \\
    --6merpair \\
    --addfragments input.fasta \\
    $ref_fasta > mafft.fasta
  """
}
