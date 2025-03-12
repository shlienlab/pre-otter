# Adding a platform tag to ensure that images built on ARM-based machines
ARG PLATFORM=linux/amd64
FROM --platform=$PLATFORM continuumio/miniconda3:23.5.2-0-alpine

ENV TERM=xterm-256color

USER root

# Install cool stuff
RUN rm -rf /var/cache/apk/* && apk upgrade && apk update & apk add git tar
    
# Update conda and Set channels for bioconda
RUN conda install conda=24.7.1 && \
    conda config --add channels nanoporetech && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda install -c bioconda star && \
    conda install -c bioconda salmon=1.10.3 && \
    conda install -y fastcat=0.20.0 && \
    conda install -y pychopper=2.7.10 && \
    conda install -y minimap2=2.28 && \
    conda install -y samtools && \
    conda install -y rsem && \
    conda install -y fastp && \
    conda install -y pip && \
    pip install requests pandas

COPY pre-otter.py /opt/
COPY ./scripts/* /opt/scripts/
COPY ./ref/* /opt/
COPY ./test_data /test_data/
RUN chmod +x /opt/pre-otter.py && ln -s /opt/pre-otter.py /usr/local/bin/pre-otter

CMD ["/bin/bash"]
