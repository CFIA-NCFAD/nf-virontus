process MULTIQC {
  publishDir "${params.outdir}/MultiQC", mode: params.publish_dir_mode

  input:
  path(multiqc_config)
  path('samtools/*')
  path('mosdepth/*')
  path('bcftools/*')
  path('software_versions/*')
  path(workflow_summary)

  output:
  path "*multiqc_report.html", emit: multiqc_report
  path "*_data"
  path "multiqc_plots"

  script:
  custom_runName = params.name
  if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
  }
  rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
  rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
  custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
  // TODO nf-core: Specify which MultiQC modules to use with -m for a faster run time
  """
  multiqc -f $rtitle $rfilename $custom_config_file .
  """
}
