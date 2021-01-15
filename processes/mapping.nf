process MAP {
  tag "$sample VS $ref_name"
  label 'process_medium'
  publishDir "${params.outdir}/mapping/$sample",
    mode: params.publish_dir_mode
  
  input:
  tuple val(sample),
        path(fastq),
        path(ref_fasta)
  output:
  tuple val(sample),
        path(ref_fasta),
        path('*.{bam,bam.bai}'), emit: bam
  path '*.{flagstat,idxstats,stats}', emit: stats
  tuple val(sample), path('*-depths.tsv'), emit: depths

  script:
  ref_name = ref_fasta.getBaseName()
  prefix = "${sample}-${ref_name}"
  bam = "${prefix}.bam"
  """
  minimap2 \\
    -ax map-ont \\
    -t${task.cpus} \\
    $ref_fasta \\
    $fastq \\
    | samtools sort -@${task.cpus} \\
    | samtools view -F4 -b -o $bam -
  samtools index $bam
  samtools stats $bam > ${prefix}.stats
  samtools flagstat $bam > ${prefix}.flagstat
  samtools idxstats $bam > ${prefix}.idxstats
  samtools depth -a -d 0 $bam | perl -ne 'chomp \$_; print "${sample}\t\$_\n"' > ${prefix}-depths.tsv
  """
}