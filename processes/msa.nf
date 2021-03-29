process MAFFT_MSA {
  label 'process_medium'

  conda (params.enable_conda ? "bioconda::mafft=7.475" : null)
  if (workflow.containerEngine == 'singularity') {
    container "https://depot.galaxyproject.org/singularity/mafft:7.475--h779adbc_1"
  } else {
    container 'quay.io/biocontainers/mafft:7.475--h779adbc_1'
  }

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
