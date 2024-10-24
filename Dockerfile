# Adding a platform tag to ensure that images built on ARM-based machines
ARG PLATFORM=linux/amd64
FROM --platform=$PLATFORM continuumio/miniconda3:23.5.2-0-alpine

ENV TERM=xterm-256color

USER root

## Install cool stuff
RUN rm -rf /var/cache/apk/* && apk upgrade && apk update & apk add git tar

# Update conda and Set channels for bioconda
RUN conda install conda=24.7.1 && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda install -c bioconda star && \
    conda install -y rsem && \
    conda install -y fastp && \
    conda install -y pip && \
    pip install requests pandas

COPY pre-otter.py /opt/
COPY ./ref/* /opt/
RUN ln -s /opt/pre-otter.py /usr/local/bin/pre-otter

CMD ["/bin/bash"]
