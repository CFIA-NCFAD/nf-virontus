// Clair3 - Symphonizing pileup and full-alignment for high-performance long-read variant calling
// GitHub: https://github.com/HKU-BAL/Clair3
// Publication: https://www.biorxiv.org/content/10.1101/2021.12.29.474431v1
process CLAIR3 {
  tag "${meta.id}"
  label 'process_low'

  conda 'bioconda::clair3==1.0.4'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/clair3:1.0.4--py39hf5e1c6e_3'
  } else {
    container 'quay.io/biocontainers/clair3:1.0.4--py39hf5e1c6e_3'
  }


  input:
  tuple val(meta), path(bam), path(bai)
  path(ref_fasta)
  // optional model_path
  path model_path

  output:
  tuple val(meta), path(vcf_fix), emit: vcf
  path(clair3_dir), emit: output_dir
  path(clair3_log), emit: log
  path("versions.yml"), emit: versions

  script:
  def prefix   = "${meta.id}"
  vcf          = "${prefix}.clair3.vcf.gz"
  vcf_fix      = "${prefix}.clair3.fix.vcf.gz"
  clair3_dir   = "${prefix}.clair3"
  clair3_log   = "${clair3_dir}/run_clair3.log"
  model_suffix = "models/${params.clair3_variant_model}"
  using_conda = (workflow.containerEngine == null || workflow.containerEngine == '')
  """
  CLAIR_BIN_DIR=\$(dirname \$(which run_clair3.sh))
  if [[ "${params.clair3_user_variant_model}" != "" ]] ; then
      MODEL_PATH=${model_path}
  else
      if [[ ${using_conda} = true ]] ; then
          MODEL_PATH="\$CLAIR_BIN_DIR/${model_suffix}"
      else [[ ${using_conda} = false ]]
          MODEL_PATH="/usr/local/bin/models/${params.clair3_variant_model}"
      fi
  fi
  
  samtools faidx $ref_fasta

  run_clair3.sh \\
    --bam_fn=$bam \\
    --ref_fn=$ref_fasta \\
    --model_path="\$MODEL_PATH"\\
    --threads=${task.cpus} \\
    --platform="ont" \\
    --output=${clair3_dir} \\
    --haploid_sensitive \\
    --enable_long_indel \\
    --fast_mode \\
    --include_all_ctgs

  CLAIR3_VERSION=\$(head -n1 ${clair3_dir}/run_clair3.log | sed 's/^.*CLAIR3 VERSION: v//; s/ .*\$//')
  zcat ${clair3_dir}/merge_output.vcf.gz \\
    | awk -v CV="\$CLAIR3_VERSION" '{ if(NR==1) { print \$0; print "##source=clair3"; print "##clair3_version=" CV; } else { print \$0; } }' \\
    | gzip -c - > ${vcf}
  # add FORMAT/SAMPLE results (AF, DP most importantly) to INFO column in VCF
  fix_clair3_vcf_info.py ${vcf} | gzip -c - > ${vcf_fix}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      clair3: \$CLAIR3_VERSION
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}
