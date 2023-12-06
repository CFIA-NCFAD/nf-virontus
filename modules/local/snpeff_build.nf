process SNPEFF_BUILD {
    tag "$fasta"
    label 'process_low'

    conda 'bioconda::snpeff=5.0'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/snpeff:5.0--0'
    } else {
        container 'quay.io/biocontainers/snpeff:5.0--0'
    }

    input:
    path fasta
    path gff

    output:
    path('snpeff_db'), emit: db
    path('*.config'), emit: config
    path('versions.yml'), emit: versions

    script:
    def basename  = fasta.baseName
    def avail_mem = 4
    if (!task.memory) {
      log.info '[snpEff] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
      avail_mem = task.memory.giga
    }
    """
    mkdir -p snpeff_db/genomes/
    cd snpeff_db/genomes/
    ln -s ../../$fasta ${basename}.fa
    cd ../../
    mkdir -p snpeff_db/${basename}/
    cd snpeff_db/${basename}/
    ln -s ../../$gff genes.gff
    cd ../../
    echo "${basename}.genome : ${basename}" > snpeff.config

    snpEff \\
        -Xmx${avail_mem}g \\
        build \\
        -config snpeff.config \\
        -dataDir ./snpeff_db \\
        -gff3 \\
        -v \\
        ${basename}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(snpEff -version 2>&1 | sed 's/^.*SnpEff //; s/ .*\$//')
    END_VERSIONS
    """
}
