process MAKE_SNPEFF_DB {
  tag "${index_base}.fa"
  label 'process_low'
  if (params.save_reference) {
      publishDir "${params.outdir}/genome", mode: params.publish_dir_mode
  }

  input:
  tuple val(index_base),
        path("SnpEffDB/genomes/${index_base}.fa"),
        path("SnpEffDB/${index_base}/genes.gff")

  output:
  tuple val(index_base),
        path("SnpEffDB"),
        path("*.config")

  script:
  """
  echo "${index_base}.genome : ${index_base}" > snpeff.config
  snpEff build -config snpeff.config -dataDir ./SnpEffDB -gff3 -v ${index_base}
  """
}

process SNPEFF {
  tag "$sample"
  label 'process_medium'
  publishDir "${params.outdir}/variants/snpeff",
    mode: params.publish_dir_mode

  input:
  tuple val(sample),
        path(ref_fasta),
        path(vcf)
  tuple val(index_base),
        path(db),
        path(config)

  output:
  path "${sample}.snpEff.csv", emit: csv
  path "*.vcf.gz*", emit: vcfgz
  path "*.{txt,html}", emit: report
  tuple val(sample), path("*.tsv"), emit: table

  script:
  """
  snpEff ${index_base} \\
    -config $config \\
    -dataDir $db \\
    ${vcf} \\
    -csvStats ${sample}.snpEff.csv \\
    | bgzip -c > ${sample}.snpEff.vcf.gz
  tabix -p vcf -f ${sample}.snpEff.vcf.gz
  mv snpEff_summary.html ${sample}.snpEff.summary.html
  SnpSift extractFields -s "," \\
    -e "." \\
    ${sample}.snpEff.vcf.gz \\
    CHROM POS REF ALT DP AC \\
    "ANN[*].GENE" "ANN[*].GENEID" \\
    "ANN[*].IMPACT" "ANN[*].EFFECT" \\
    "EFF[*].CODON" "EFF[*].AA" \\
    "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \\
    "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \\
    "ANN[*].AA_LEN" \\
    > ${sample}.snpSift.table.txt
  simpler_snpsift.py ${sample}.snpSift.table.txt ${sample}.snpSift.simple.tsv
  """
}
