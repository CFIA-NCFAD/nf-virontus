process SAMTOOLS_SORT_FGBIO {
  tag "$meta.id"
  label "process_low"

  conda 'bioconda::samtools=1.16.1'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1'
  } else {
    container 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'
  }

  input:
  tuple val(meta), path(bam), path(bai)
  

  output:
  tuple val(meta), path('*.sort.bam'), emit: bam
  path('versions.yml'), emit: versions

  script:
  def prefix = "${meta.id}.sort"
  """
  samtools sort -n -u $bam > ${prefix}.bam

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}
