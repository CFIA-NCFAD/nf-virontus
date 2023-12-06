process IVAR_TRIM {
  tag "$meta.id"
  label "process_low"

  conda 'bioconda::ivar=1.4.2'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/ivar:1.4.2--h0033a41_2'
  } else {
    container 'quay.io/biocontainers/ivar:1.4.2--h0033a41_2'
  }

  input:
  tuple val(meta), path(bam), path(bai)
  path(bed)

  output:
  tuple val(meta), path('*.trim.bam'), path('*.bam.bai'), emit: bam
  path('*.{flagstat,idxstats,stats}'),                    emit: stats
  tuple val(meta), path('*.samtools.depth'),              emit: depth
  path('versions.yml'),                                   emit: versions

  script:
  def prefix = "${meta.id}.trim"
  def keep_reads_no_primer = params.ivar_trim_noprimer ? '' : '-e'
  """
  ivar trim \\
    -i $bam \\
    -b $bed \\
    $keep_reads_no_primer \\
    -q 1 -m 20 -s 4 \\
    -p trim
  samtools sort -o ${prefix}.bam trim.bam
  rm trim.bam
  samtools index ${prefix}.bam
  samtools stats ${prefix}.bam > ${prefix}.stats
  samtools flagstat ${prefix}.bam > ${prefix}.flagstat
  samtools idxstats ${prefix}.bam > ${prefix}.idxstats
  samtools depth -a -d 0 ${prefix}.bam > ${prefix}.samtools.depth

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      ivar: \$(ivar version | head -n1 | sed 's/iVar version //')
  END_VERSIONS
  """
}
