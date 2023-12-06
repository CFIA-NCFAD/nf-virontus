process NEXTCLADE_DATASETGET {
  label 'process_low'

  conda "bioconda::nextclade=2.14.0"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/nextclade:2.14.0--h9ee0642_1'
  } else {
    container 'quay.io/biocontainers/nextclade:2.14.0--h9ee0642_1'
  }

  input:
  val(dataset)

  output:
  path("nextclade-dataset/"), emit: nextclade_dataset
  path("versions.yml"), emit: versions

  script:
  """
  nextclade dataset get \\
    --name "${dataset}" \\
    --output-dir nextclade-dataset

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      nextclade: \$(nextclade --version)
  END_VERSIONS
  """
}


process NEXTCLADE_RUN {
  label 'process_high_cpu_medium_mem'

  conda "bioconda::nextclade=2.14.0"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/nextclade:2.14.0--h9ee0642_1'
  } else {
    container 'quay.io/biocontainers/nextclade:2.14.0--h9ee0642_1'
  }

  input:
  path(fasta)
  val(dataset)

  output:
  path("nextclade.tsv")                , optional:true, emit: tsv
  path("nextclade.csv")                , optional:true, emit: csv
  path("nextclade.json")               , optional:true, emit: json
  path("nextclade.auspice.json")       , optional:true, emit: auspice_json
  path("nextclade-output/")            , optional:true, emit: nextclade_output
  path("versions.yml")                 , emit: versions

  script:
  """
  nextclade run \\
    --jobs $task.cpus \\
    --output-all nextclade-output \\
    --input-dataset $dataset \\
    --output-tsv nextclade.tsv \\
    --output-csv nextclade.csv \\
    --output-tree nextclade.auspice.json \\
    --output-basename nextclade \\
    $fasta

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      nextclade: \$(nextclade --version)
  END_VERSIONS
  """
}
