FROM mambaorg/micromamba:2.0.5-debian12-slim

# Set up environment variables
ENV PATH=/usr/local/bin:$PATH
USER root

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get -y install gcc mono-mcs && \
    rm -rf /var/lib/apt/lists/*

# Install GATK, Samtools, BWA, and Tabix (htslib)
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.10 \
    "gatk4>=4.2.6" \
    "bwa>=0.7.17" \
    samtools \
    htslib \
    biopython \
    pysam \
    numpy \
    matplotlib \
    scipy \
    pandas \
    h5py \
    && micromamba clean -a -y

# Test pysam installation
RUN micromamba run -n base python3 -c "import pysam; print(pysam.__version__)"

# Copy and install DupCaller
COPY setup.py /opt/DupCaller/setup.py
COPY src /opt/DupCaller/src
WORKDIR /opt/DupCaller
RUN micromamba run -n base pip install -e .

CMD ["bash"]