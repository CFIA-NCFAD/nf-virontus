process BCFTOOLS_STATS {
  tag "${meta.id}"
  
  conda "bioconda::bcftools=1.18"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/bcftools:1.18--h8b25389_0'
  } else {
    container 'quay.io/biocontainers/bcftools:1.18--h8b25389_0'
  }

  input:
  tuple val(meta), path(vcf)
  path(ref_fasta)

  output:
  path("*.bcftools_stats.txt"), emit: stats
  path('versions.yml'), emit: versions

  script:
  """
  bcftools stats -F $ref_fasta $vcf > ${meta.id}.bcftools_stats.txt

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
  END_VERSIONS
  """
}
