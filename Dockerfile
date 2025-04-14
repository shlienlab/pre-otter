# Adding a platform tag to ensure that images built on ARM-based machines
ARG PLATFORM=linux/amd64
FROM --platform=$PLATFORM continuumio/miniconda3:23.5.2-0-alpine

ENV TERM=xterm-256color

USER root

# Install cool stuff
RUN rm -rf /var/cache/apk/* && apk upgrade && apk update && apk add \
    git \
    tar \
    libxml2 \
    openjdk11 \
    build-base \
    wget \
    g++ \
    gcc \
    make \
    curl \
    R \
    coreutils

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

# Clone JAFFA repository and checkout the required commit
RUN git clone https://github.com/Oshlack/JAFFA.git /JAFFA && \
    cd /JAFFA && \
    git checkout 24b1c3b

# Install tools and set up environment
RUN mkdir -p /JAFFA/tools/bin && cd /JAFFA/tools && \
    commands="bpipe reformat extract_seq_from_fasta make_simple_read_table process_transcriptome_align_table make_3_gene_fusion_table dedupe" && \
    # bpipe installation
    wget -O bpipe-0.9.9.2.tar.gz https://github.com/ssadedin/bpipe/releases/download/0.9.9.2/bpipe-0.9.9.2.tar.gz && \
    tar -zxvf bpipe-0.9.9.2.tar.gz && \
    rm bpipe-0.9.9.2.tar.gz && \
    ln -s /JAFFA/tools/bpipe-0.9.9.2/bin/* /JAFFA/tools/bin/ && \
    # Custom installs for each tool
    g++ -std=c++11 -O3 -o /JAFFA/tools/bin/make_3_gene_fusion_table /JAFFA/src/make_3_gene_fusion_table.c++ && \
    g++ -std=c++11 -O3 -o /JAFFA/tools/bin/extract_seq_from_fasta /JAFFA/src/extract_seq_from_fasta.c++ && \
    g++ -std=c++11 -O3 -o /JAFFA/tools/bin/make_simple_read_table /JAFFA/src/make_simple_read_table.c++ && \
    g++ -std=c++11 -O3 -o /JAFFA/tools/bin/process_transcriptome_align_table /JAFFA/src/process_transcriptome_align_table.c++ && \
    g++ -O3 -o /JAFFA/tools/bin/make_count_table /JAFFA/src/make_count_table.c++ && \
    # Dedupe installation
    wget --no-check-certificate https://sourceforge.net/projects/bbmap/files/BBMap_36.59.tar.gz && \
    tar -zxvf BBMap_36.59.tar.gz && \
    rm BBMap_36.59.tar.gz && \
    for script in $(ls /JAFFA/tools/bbmap/*.sh); do \
        s=$(basename $script); \
        s_pre=$(echo $s | sed 's/.sh//g'); \
        echo "/JAFFA/tools/bbmap/$s \$@" > /JAFFA/tools/bin/$s_pre; \
        chmod +x /JAFFA/tools/bin/$s_pre; \
    done

# Check if GCC version is >= 4.9
RUN gcc_version=$(gcc -dumpversion) && \
    gcc_check=$(echo -e "$gcc_version\n4.9" | sort -V | tail -n1) && \
    if [ "$gcc_check" != "$gcc_version" ]; then \
        echo "gcc must be >= 4.9. Exiting..."; \
        exit 1; \
    fi

# Verify the tools directory contents
RUN ls /JAFFA/tools/bin/

# Final check to ensure all required tools were installed
RUN echo "Checking that all required tools were installed:" && \
    Final_message="All commands installed successfully!" && \
    for c in $commands; do \
        c_path=$(which /JAFFA/tools/bin/$c 2>/dev/null); \
        if [ -z "$c_path" ]; then \
            echo -n "WARNING: $c could not be found!!!! "; \
            echo "You will need to download and install $c manually."; \
            Final_message="WARNING: One or more commands did not install successfully."; \
        else \
            echo "$c looks like it has been installed"; \
        fi; \
    done && \
    echo "**********************************************************" && \
    echo "$Final_message"

RUN cp /opt/conda/bin/minimap2 /JAFFA/tools/bin/
ENV PATH=/JAFFA/tools/bin:$PATH

# Create tools.groovy and populate it with the correct tool paths
RUN echo "// Path to tools used by the JAFFA pipeline" > /JAFFA/tools.groovy && \
    commands="bpipe reformat extract_seq_from_fasta make_simple_read_table process_transcriptome_align_table make_3_gene_fusion_table dedupe" && \
    for c in $commands; do \
        c_path=$(which /JAFFA/tools/bin/$c 2>/dev/null); \
        if [ -z "$c_path" ]; then \
            echo "$c not found, fetching it"; \
            /JAFFA/tools/$c"_install"; \
            c_path=$(which /JAFFA/tools/bin/$c 2>/dev/null); \
        fi; \
        echo "$c=\"$c_path\"" >> /JAFFA/tools.groovy; \
    done && \
    # Check if R is installed
    R_path=$(which R 2>/dev/null) && \
    if [ -z "$R_path" ]; then \
        echo "R not found! Please install it manually."; \
    fi && \
    echo "R=\"$R_path\"" >> /JAFFA/tools.groovy && \
    # Check if minimap2 is installed
    minimap2_path=$(which minimap2 2>/dev/null) && \
    if [ -z "$minimap2_path" ]; then \
        echo "minimap2 not found! Please install it manually."; \
    else \
        echo "minimap2=\"$minimap2_path\"" >> /JAFFA/tools.groovy; \
    fi

COPY pre-otter.py /opt/
COPY ./scripts/* /opt/scripts/
COPY ./ref/* /opt/
COPY ./test_data /test_data/
RUN chmod +x /opt/scripts/download_jaffal_references.sh
RUN chmod +x /opt/pre-otter.py && ln -s /opt/pre-otter.py /usr/local/bin/pre-otter

CMD ["/bin/bash"]
