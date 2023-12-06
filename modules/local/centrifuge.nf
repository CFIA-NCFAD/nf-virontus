process CENTRIFUGE {
  tag "$meta.id"
  label 'process_highmem'

  conda "bioconda::centrifuge-core=1.0.4"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/centrifuge-core:1.0.4--h43eeafb_2'
  } else {
    container 'quay.io/biocontainers/centrifuge-core:1.0.4--h43eeafb_2'
  }

  input:
  tuple path(db), val(db_prefix)
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path("*.centrifuge.report.txt")  , emit: report
  tuple val(meta), path("*.centrifuge.results.txt") , emit: results
  path  "*.version.txt"                             , emit: version

  script:
  def software = getSoftwareName(task.process)
  def prefix = "${meta.id}"
  """
  centrifuge \\
    -x ${db}/${db_prefix} \\
    -U $reads \\
    -S ${prefix}.centrifuge.results.txt \\
    ${option.args} \\
    -p ${task.cpus}

  centrifuge-kreport \\
    -x ${db}/${db_prefix} \\
    ${prefix}.centrifuge.results.txt \\
    > ${prefix}.centrifuge.report.txt

  centrifuge --version | head -n1 | sed -E 's/^.*centrifuge-class version (\\S+)/\\1/' > ${software}.version.txt
  """
}
