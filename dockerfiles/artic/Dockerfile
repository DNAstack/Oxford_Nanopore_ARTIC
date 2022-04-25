# TAG dnastack/artic:8af5152
FROM continuumio/miniconda3:4.10.3 AS condabuild

MAINTAINER Heather Ward <heather@dnastack.com>

RUN git clone https://github.com/connor-lab/ncov2019-artic-nf.git /opt/ncov2019-artic-nf

WORKDIR /opt/ncov2019-artic-nf

ENV ARTIC_VERSION 8af5152cf7107c3a369b996f5bad3473d329050c
RUN git checkout ${ARTIC_VERSION}

RUN /opt/conda/bin/conda install mamba -c conda-forge && \
    /opt/conda/bin/mamba env create -f environments/nanopore/environment.yml

RUN apt-get -qq update && \
    apt-get -qq install \
        tabix \
        samtools \
        xxd \
        zip

COPY requirements.txt /opt/requirements.txt
RUN python3 -m pip install -r /opt/requirements.txt

ENV SCHEME_REPO /opt/artic-ncov2019

RUN git clone https://github.com/artic-network/artic-ncov2019.git "$SCHEME_REPO" && \
    for reference in $(find "${SCHEME_REPO}/primer_schemes/nCoV-2019" -name "nCoV-2019.reference.fasta"); do \
        samtools faidx "${reference}"; \
    done

COPY ./scripts/* /usr/local/bin/

ENV PATH=/opt/conda/envs/artic/bin:/opt/ncov2019-artic-nf/bin:$PATH
