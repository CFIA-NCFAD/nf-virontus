process MOSDEPTH_GENOME {
  tag "$sample"
  label 'process_medium'
  publishDir "${params.outdir}/mosdepth", 
             mode: params.publish_dir_mode,
             saveAs: { filename ->
               if (filename.endsWith(".pdf")) "plots/$filename"
               else if (filename.endsWith(".tsv")) "plots/$filename"
               else filename
             }

  input:
  tuple val(sample), path(bam)

  output:
  path "*.global.dist.txt", emit: mqc
  path "*.{txt,gz,csi,tsv}"

  script:
  """
  mosdepth \\
      --fast-mode \\
      $sample \\
      ${bam[0]}
  """
}
