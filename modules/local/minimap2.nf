process MINIMAP2 {
  tag "${meta.id}"
  label 'process_medium'

  conda 'bioconda::minimap2=2.18 bioconda::samtools=1.12'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:836eb07132d2de763323bbd4f1083b3fdf328759-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:836eb07132d2de763323bbd4f1083b3fdf328759-0'
  }

  input:
  tuple val(meta), path(fastq), path(ref_fasta)

  output:
  tuple val(meta), path('*.bam'), path('*.bam.bai'), emit: bam
  path('*.{flagstat,idxstats,stats}'),               emit: stats
  tuple val(meta), path('*.samtools.depth'),         emit: depth
  path('versions.yml'),                              emit: versions

  script:
  def ref_name = ref_fasta.baseName
  def prefix = "${meta.id}"
  bam = "${meta.id}.bam"
  """
  minimap2 \\
    -ax map-ont \\
    -t${task.cpus} \\
    $ref_fasta \\
    $fastq \\
    | samtools sort -@${task.cpus} \\
    > $bam

  samtools index $bam
  samtools stats $bam > ${prefix}.stats
  samtools flagstat $bam > ${prefix}.flagstat
  samtools idxstats $bam > ${prefix}.idxstats
  samtools depth -a -d 0 ${prefix}.bam > ${prefix}.samtools.depth

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      minimap2: \$(minimap2 --version)
      samtools: \$(samtools --version 2>&1 | tr '\\n' '\\t' | sed  's/.*samtools //; s/\\t.*//')
  END_VERSIONS
  """
}
