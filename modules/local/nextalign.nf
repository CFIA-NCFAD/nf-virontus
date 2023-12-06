process NEXTALIGN {
  label 'process_high'

  conda 'bioconda::nextalign=2.7.0'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/nextalign:2.7.0--h9ee0642_0'
  } else {
    container 'quay.io/biocontainers/nextalign:2.7.0--h9ee0642_0'
  }

  input:
  path(fastas)
  path(ref_fasta)

  output:
  path('nextalign.fasta'),          emit: msa
  path('nextalign.insertions.csv'), emit: insertions
  path('versions.yml'),             emit: versions

  script:
  """
  # keep only the seq id; remove all other text after fasta seq id
  cat $fastas | sed -r 's/^>(\\S+)\\s.*/>\\1/g' > input.fasta
  nextalign \\
    -j ${task.cpus} \\
    --sequences=input.fasta \\
    --reference=$ref_fasta \\
    --output-dir=./ \\
    --output-basename=nextalign

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      nextalign: \$(nextalign --version | sed 's/nextalign //')
  END_VERSIONS
  """
}
