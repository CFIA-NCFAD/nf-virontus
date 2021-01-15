process IVAR_TRIM {
  tag "$sample - $ref_name"
  label "process_low"

  publishDir "${params.outdir}/mapping/$sample",
    mode: params.publish_dir_mode

  input:
  path(bedfile)
  tuple val(sample),
        path(ref_fasta),
        path(bam)
  output:
  tuple val(sample),
        path(ref_fasta),
        path('*.trim.{bam,bam.bai}'), emit: bam
  path '*.{flagstat,idxstats,stats}', emit: stats
  tuple val(sample), path('*-depths.tsv'), emit: depths

  script:
  ref_name = ref_fasta.getBaseName()
  prefix = "${sample}-${ref_name}.trim"
  """
  ivar trim \\
    -i ${bam[0]} \\
    -b $bedfile \\
    -p trim -q 1 -m 20 -s 4 -e
  samtools sort -o ${prefix}.bam trim.bam
  samtools index ${prefix}.bam
  rm trim.bam
  samtools stats ${prefix}.bam > ${prefix}.stats
  samtools flagstat ${prefix}.bam > ${prefix}.flagstat
  samtools idxstats ${prefix}.bam > ${prefix}.idxstats
  samtools depth -a -d 0 ${prefix}.bam | perl -ne 'chomp \$_; print "${sample}\t\$_\n"' > ${prefix}-depths.tsv
  """
}