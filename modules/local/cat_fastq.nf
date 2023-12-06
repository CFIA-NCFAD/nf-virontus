// Concatenate FASTQ files from the same sample into single file and gzip compress with pigz if not already gzipped.
process CAT_FASTQ {
  tag "${meta.id}"
  label 'process_medium'

  conda 'conda-forge::pigz=2.6'
  // use smallest pre-built container with pigz 2.6 and prodigal (7MB)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-2e442ba7b07bfa102b9cf8fac6221263cd746ab8:57f05cfa73f769d6ed6d54144cb3aa2a6a6b17e0-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-2e442ba7b07bfa102b9cf8fac6221263cd746ab8:57f05cfa73f769d6ed6d54144cb3aa2a6a6b17e0-0'
  }

  input:
  tuple val(meta), path(fqgz), path(fq)

  output:
  tuple val(meta), path(merged_fqgz), emit: reads
  path "versions.yml", emit: versions

  script:
  merged_fqgz = "${meta.id}.merged.fastq.gz"
  def fqList = fq.collect { it.toString() }
  def fqgzList = fqgz.collect { it.toString() }
  """
  touch $merged_fqgz
  if [ ${fqList.size} -gt 0 ]; then
    cat $fq | pigz -ck >> $merged_fqgz
  fi
  if [ ${fqgzList.size} -gt 0 ]; then
    cat $fqgz >> $merged_fqgz
  fi
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      pigz: \$(pigz --version 2>&1 | sed 's/pigz //g' )
  END_VERSIONS
  """
}
