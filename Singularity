Bootstrap:docker
From:continuumio/miniconda3:4.7.12

%labels
    MAINTAINER Peter Kruczkiewicz
    DESCRIPTION Singularity image containing all requirements for the peterk87/nf-virontus pipeline
    VERSION 1.0.0

%environment
    PATH=/opt/conda/envs/nf-virontus-1.0.0/bin:$PATH
    PATH=/opt/guppy/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    export PATH=/opt/conda/bin:$PATH
    apt-get update && apt-get install -y procps curl && apt-get clean -y
    curl -SsL -o guppy.tgz https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_3.4.4_linux64.tar.gz
    mkdir /opt/guppy
    tar -C /opt/guppy --strip-components=1 -xf guppy.tgz
    rm guppy.tgz
    conda update conda
    conda env create -f /environment.yml
    conda clean -a
