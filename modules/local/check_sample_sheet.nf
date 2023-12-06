process CHECK_SAMPLE_SHEET {
  conda 'bioconda::shiptv=0.4.1'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/shiptv:0.4.1--pyh5e36f6f_0'
  } else {
    container 'quay.io/biocontainers/shiptv:0.4.1--pyh5e36f6f_0'
  }

  input:
  path samplesheet

  output:
  path "samplesheet_reformat.csv" 

  script:
  """
  check_sample_sheet.py $samplesheet samplesheet_reformat.csv
  """
}
