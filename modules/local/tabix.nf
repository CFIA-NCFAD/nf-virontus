process TABIX {
  tag "$meta.id"
  
  conda "bioconda::tabix=0.2.6"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/tabix:0.2.6--ha92aebf_0"
  } else {
    container "quay.io/biocontainers/tabix:0.2.6--ha92aebf_0"
  }

  input:
  tuple val(meta), path(vcf)

  output:
  tuple val(meta), path('*.vcf.gz'), path('*.tbi'), emit: tbi
  path('versions.yml'), emit: versions

  script:
  """
  bgzip -c $vcf > ${vcf.name}.gz
  tabix -p vcf -f ${vcf.name}.gz

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      tabix: \$(tabix -h 2>&1 | tr '\n' '\t' | sed 's/.*Version: //; s/\t.*//')
  END_VERSIONS
  """
}
