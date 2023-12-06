process COVERAGE_PLOT {
  tag "$meta.id"
  publishDir "${params.outdir}/plots/coverage",
    mode: params.publish_dir_mode

  conda 'python=3.9 conda-forge::typer=0.3.2 conda-forge::rich=10.6.0 conda-forge::seaborn=0.11.0 conda-forge::pandas=1.3.0 bioconda::bcbio-gff=0.6.6 bioconda::dna_features_viewer=3.0.3'

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-596f42d854e849eb773ecd1b48f2b698c2d09c9f:400d0a2593841aa0bfa3402fe85debd55a29cf37-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-596f42d854e849eb773ecd1b48f2b698c2d09c9f:400d0a2593841aa0bfa3402fe85debd55a29cf37-0'
  }


  input:
  tuple val(meta), path(vcf), path(depths)
  path(ref_fasta)
  
  output:
  path("*.pdf")

  script:
  def ref_name = ref_fasta.getBaseName()
  def prefix = "${meta.id}.${ref_name}.coverage"
  def sample = meta.id
  """
  # coverage plot with variants
  plot_coverage.py \\
    -d $depths \\
    -v $vcf \\
    --sample-name $sample \\
    --no-highlight \\
    --low-coverage ${params.low_coverage} \\
    -o ${prefix}.variants.pdf

  plot_coverage.py \\
    -d $depths \\
    -v $vcf \\
    --sample-name $sample \\
    --log-scale-y \\
    --no-highlight \\
    --low-coverage ${params.low_coverage} \\
    -o ${prefix}.variants.log-scale-y.pdf

  # coverage plot highlighting low/no coverage regions
  plot_coverage.py \\
    -d $depths \\
    --sample-name $sample \\
    --log-scale-y \\
    --low-coverage ${params.low_coverage} \\
    -o ${prefix}.log-scale-y.highlight-no-low-cov.pdf

  plot_coverage.py \\
    -d $depths \\
    --sample-name $sample \\
    --low-coverage ${params.low_coverage} \\
    -o ${prefix}.highlight-no-low-cov.pdf

  # coverage plot no highlighting
  plot_coverage.py \\
    -d $depths \\
    --sample-name $sample \\
    --log-scale-y \\
    --no-highlight \\
    --low-coverage ${params.low_coverage} \\
    -o ${prefix}.log-scale-y.pdf

  plot_coverage.py \\
    -d $depths \\
    --sample-name $sample \\
    --no-highlight \\
    --low-coverage ${params.low_coverage} \\
    -o ${prefix}.pdf
  """
}

process BASIC_TREE_PLOT {
  publishDir "${params.outdir}/plots/tree",
             mode: params.publish_dir_mode
  input:
  path(iqtree_output)
  path(pangolin_report)

  output:
  path(output)

  script:
  output = 'IQ-TREE_Maximum-Likelihood_Phylogenetic_Tree_mqc.png'
  pangolin_report_arg = pangolin_report ? "--pangolin-report $pangolin_report" : ""
  """
  plot_basic_tree.py tree-*.treefile $output $pangolin_report_arg
  """
}
