process SNPSIFT {
  tag "$meta.id"
  label 'process_low'

  conda 'bioconda::snpsift=4.3.1t'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/snpsift:4.3.1t--hdfd78af_3'
  } else {
    container 'quay.io/biocontainers/snpsift:4.3.1t--hdfd78af_3'
  }

  input:
  tuple val(meta), path(snpeff_vcf_gz), path(tabix_tbi)

  output:
  tuple val(meta), path('*.snpsift.table.txt'), emit: txt
  path('versions.yml'), emit: versions

  script:
  """
  SnpSift extractFields -s "," \\
    -e "." \\
    $snpeff_vcf_gz \\
    CHROM POS REF ALT DP AF \\
    "ANN[*].GENE" "ANN[*].GENEID" \\
    "ANN[*].IMPACT" "ANN[*].EFFECT" \\
    "EFF[*].CODON" "EFF[*].AA" \\
    "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \\
    "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \\
    "ANN[*].AA_LEN" \\
    > ${meta.id}.snpsift.table.txt

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      snpsift: \$(SnpSift sort -h 2>&1 | tr '\\n' '\\t' | sed 's/.*SnpSift version //;s/ .build .*//')
  END_VERSIONS
  """
}

