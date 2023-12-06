process VCF_FILTER_FRAMESHIFT {
  tag "$meta.id"
  label 'process_low'
  
  conda 'bioconda::shiptv=0.4.1'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/shiptv:0.4.1--pyh5e36f6f_0'
  } else {
    container 'quay.io/biocontainers/shiptv:0.4.1--pyh5e36f6f_0'
  }

  input:
  tuple val(meta), path(vcf)

  output:
  tuple val(meta), path(out_vcf)

  script:
  out_vcf = "${file(vcf).baseName}.no_fs.vcf"
  """
  vcf_filter_frameshift.py $vcf $out_vcf
  """
}
