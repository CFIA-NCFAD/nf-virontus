process SIMPLER_SNPSIFT {
  tag "$meta.id"

  conda 'bioconda::shiptv=0.4.1'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/shiptv:0.4.1--pyh5e36f6f_0'
  } else {
    container 'quay.io/biocontainers/shiptv:0.4.1--pyh5e36f6f_0'
  }

  input:
  tuple val(meta), path(snpsift_table)

  output:
  tuple val(meta), path('*.snpsift.simple.tsv'), emit: tsv

  script:
  """
  simpler_snpsift.py $snpsift_table ${meta.id}.snpsift.simple.tsv
  """
}
