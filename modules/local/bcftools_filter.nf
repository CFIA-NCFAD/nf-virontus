// Exclude variants that have less than the required allele fraction (AF)
// calculated from unambiguous variant calls and with less unambiguous depth
// than the low coverage threshold.
// The SR INFO tag produced by Medaka annotation is essential for this
// calculation.
process BCFTOOLS_FILTER {
  tag "$meta.id - AF$allele_fraction"
  label 'process_low'

  conda "bioconda::bcftools=1.18"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/bcftools:1.18--h8b25389_0'
  } else {
    container 'quay.io/biocontainers/bcftools:1.18--h8b25389_0'
  }

  input:
  tuple val(meta), path(vcf)
  val allele_fraction

  output:
  tuple val(meta), path(filt_vcf), emit: vcf
  path('versions.yml'), emit: versions

  script:
  filt_vcf = "${meta.id}.${allele_fraction}AF.filt.vcf"
  """
  bcftools filter \\
    -e '(AF < $allele_fraction) || (DP < ${params.low_coverage})' \\
    $vcf \\
    -Ov \\
    -o $filt_vcf

  bcftools --version | head -n1 | sed 's/^bcftools //' > bcftools.version.txt
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
  END_VERSIONS
  """
}
