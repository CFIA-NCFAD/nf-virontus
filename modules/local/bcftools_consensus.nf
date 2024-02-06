process BCFTOOLS_CONSENSUS {
  tag "$meta.id"
  label 'process_low'

  conda "bioconda::bcftools=1.18"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/bcftools:1.18--h8b25389_0'
  } else {
    container 'quay.io/biocontainers/bcftools:1.18--h8b25389_0'
  }

  input:
  tuple val(meta), path(vcf_gz), path(tabix_tbi), path(mosdepth_per_base)
  path(ref_fasta)
  val(low_coverage)

  output:
  tuple val(meta), path(consensus), emit: consensus
  path('versions.yml'), emit: versions

  script:
  consensus = "${meta.id}.consensus.fasta"
  """
  zcat $mosdepth_per_base | awk '\$4<${low_coverage}' > low_cov.bed
  bcftools consensus \\
    -H A \\
    -f $ref_fasta \\
    -m low_cov.bed \\
    $vcf_gz > $consensus

  sed -i -E "s/^>(.*)/>${meta.id}/g" $consensus

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
  END_VERSIONS
  """
}
