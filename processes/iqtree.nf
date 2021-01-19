process IQTREE {
  publishDir "${params.outdir}/iqtree",
             mode: params.publish_dir_mode

  input:
  path("mafft.fasta")
  val(model)

  output:
  path('tree-*')

  script:
  """
  OUTGROUP=`grep '^>' mafft.fasta | head -n1 | sed -r 's/>(\\S+)/\\1/'`
  iqtree \\
    -s mafft.fasta \\
    -o \$OUTGROUP \\
    -T ${task.cpus} \\
    -m $model \\
    --prefix tree-\$OUTGROUP-${model}_model
  """
}
