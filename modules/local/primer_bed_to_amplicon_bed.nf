process PRIMER_BED_TO_AMPLICON_BED {
  label 'process_low'

  conda "bioconda::shiptv=0.4.1"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/shiptv:0.4.1--pyh5e36f6f_0'
  } else {
    container 'quay.io/biocontainers/shiptv:0.4.1--pyh5e36f6f_0'
  }

  input:
  path(primer_bed)

  output:
  path(amplicon_bed)

  script:
  amplicon_bed = "${file(primer_bed).baseName}.amplicon.bed"
  """
  primer_bed_to_amplicon_bed.py \\
    $primer_bed \\
    $amplicon_bed \\
    --left-suffix ${params.left_primer_suffix} \\
    --right-suffix ${params.right_primer_suffix}
  """
}
