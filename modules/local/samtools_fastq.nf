process SAMTOOLS_FASTQ {
  tag "$meta.id"
  label 'process_low'

  conda "bioconda::samtools=1.16.1 conda-forge::pigz=2.6"
  // mulled container with bowtie2=2.4.5,samtools=1.15,pigz=2.6 but only 95MB
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:67d2f6bfef729f65a6870792b4467bd70e26b078-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:67d2f6bfef729f65a6870792b4467bd70e26b078-0'
  }

  input:
  tuple val(meta), path(bam)

  output:
  tuple val(meta), path('*.fastq.gz'), emit: reads
  path('versions.yml'), emit: versions

  script:
  def prefix = "${meta.id}.fgbio.clipbam"
  """
  samtools fastq $bam | pigz -ck > ${prefix}.fastq.gz

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
      pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
  END_VERSIONS
  """
}
