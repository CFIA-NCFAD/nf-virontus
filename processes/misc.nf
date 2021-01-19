process CHECK_SAMPLE_SHEET {
  tag "$samplesheet"
  publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

  input:
  path samplesheet

  output:
  path "samplesheet_reformat.csv" 

  script:
  """
  check_sample_sheet.py $samplesheet samplesheet_reformat.csv
  """
}
