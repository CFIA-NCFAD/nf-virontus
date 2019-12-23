FROM continuumio/miniconda3:4.7.12
LABEL authors="Peter Kruczkiewicz" \
      description="Docker image containing all requirements for peterk87/nf-virontus pipeline"

COPY environment.yml /
RUN apt-get update && \
    apt-get install -y procps curl && \
    apt-get clean -y && \
    curl -SsL -o guppy.tgz https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_3.4.4_linux64.tar.gz && \
    mkdir /opt/guppy && \
    tar -C /opt/guppy --strip-components=1 -xf guppy.tgz && \
    rm guppy.tgz && \
    conda update conda && \
    conda env create -f /environment.yml && \
    conda clean -a
ENV PATH /opt/conda/envs/nf-virontus-1.0.0/bin:$PATH
ENV PATH /opt/guppy/bin:$PATH
