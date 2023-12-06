process FGBIO_CLIPBAM {
  tag "${meta.id}"
  label 'process_low'

  conda 'bioconda::fgbio=2.0.2'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/fgbio:2.0.2--hdfd78af_0'
  } else {
    container 'quay.io/biocontainers/fgbio:2.0.2--hdfd78af_0'
  }

  input:
  tuple val(meta), path(bam)
  path(ref_fasta)

  output:
  tuple val(meta), path("*.fgbio.clipbam.bam"), emit: bam
  path('versions.yml'), emit: versions

  script:
  def prefix = "${meta.id}.fgbio.clipbam"
  """
  fgbio ClipBam \\
    -H -c Hard \\
    -r $ref_fasta \\
    -i $bam \\
    -o ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$(echo \$(fgbio --version 2>&1) | sed 's/Version: //')
    END_VERSIONS
  """
}
