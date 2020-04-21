FROM continuumio/miniconda3:4.8.2
LABEL authors="Peter Kruczkiewicz" \
      description="Docker image containing all requirements for peterk87/nf-virontus pipeline"

COPY environment.yml /
RUN apt-get update && \
    apt-get install -y procps curl && \
    apt-get clean -y && \
    conda update conda && \
    conda env create -f /environment.yml && \
    conda clean -a
ENV PATH /opt/conda/envs/nf-virontus-1.0.0/bin:$PATH
