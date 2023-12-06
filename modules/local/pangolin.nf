process PREPARE_FASTA_FOR_PANGOLIN {
  input:
  tuple val(meta), path(fasta)

  output:
  path(output)

  script:
  def sample = "${meta.id}"
  output = "${sample}.fasta"
  """
  cat $fasta > $output
  sed -i 's/^>.*/>${sample}/g' $output
  """
}

process PANGOLIN {
  label 'process_low'

  conda 'bioconda::pangolin=4.3'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/pangolin:4.3--pyhdfd78af_2'
  } else {
    container 'quay.io/biocontainers/pangolin:4.3--pyhdfd78af_2'
  }

  input:
  path(fasta)

  output:
  path 'pangolin.csv' , emit: report
  path 'pangolin.log' , emit: log
  path 'versions.yml', emit: versions

  script:
  """
  pangolin \\
    $fasta \\
    --outfile pangolin.csv \\
    --threads $task.cpus
  ln -s .command.log pangolin.log

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      pangolin: \$(pangolin --version | sed "s/pangolin //g")
  END_VERSIONS
  """
}

process PANGOLIN_SUMMARY_FOR_MULTIQC {
  conda 'bioconda::shiptv=0.4.1'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/shiptv:0.4.1--pyh5e36f6f_0'
  } else {
    container 'quay.io/biocontainers/shiptv:0.4.1--pyh5e36f6f_0'
  }

  input:
  path(input)

  output:
  path("pangolin_mqc.txt")

  script:
  """
  summarize_pangolin_results.py ./ pangolin_mqc.txt
  """
}
