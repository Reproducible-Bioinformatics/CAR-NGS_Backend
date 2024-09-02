# Usa una base di Ubuntu
FROM ubuntu:20.04

# Imposta variabili di ambiente per evitare problemi durante l'installazione
ENV DEBIAN_FRONTEND=noninteractive

# Installa i pacchetti necessari
RUN apt-get update && \
    apt-get install -y \
    curl \
    bzip2 \
    build-essential \
    zlib1g-dev \
    wget \
    python3 \
    python3-pip \
    git \
    unzip \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libhdf5-dev \
    libfftw3-dev \
    libgsl0-dev \
    libpng-dev \
    libjpeg-dev \
    libreadline-dev \
    r-base \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Installazione di Cell Ranger
RUN wget -O cellranger-8.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-8.0.1.tar.gz?Expires=1723948839&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=O~jX8XOJ3MIiqeARv0jjVX~GTvWz601rGcYBNQOw68FoREH~ZXbrNvvu4kxPTbex0NzfPwEMKs2FRaoFoQz50gMNAYu4EOvUdPUQ1Kbx10HA0hPou0h0IoTPNCp5KAj6mHHt31icMQwrpRMqz11fDh6MSKX4H8ky2q4Pg-qdFY04RrjAlHl9LhPMTit8j3gQ7FpWAh9hfvj2VL3x8tjxdTtYn4pyx-ithWgaR-kBPPDev6SM9P-b0ZDPq2NvwVUuZTP4VD5aSFsq6iP66LKc-WgvKdhlcgM14c55duNACNSx5a7cdwzjFRMSTCgjllnDqTXpNCPrjdQ5RqqT21AMlw__" && \
    tar -xzvf cellranger-8.0.1.tar.gz -C /opt/ && \
    ln -s /opt/cellranger-8.0.1/cellranger /usr/local/bin/cellranger && \
    rm cellranger-8.0.1.tar.gz

# Installazione di librerie R per l'analisi dei dati e la visualizzazione
RUN R -e "install.packages(c('ggplot2', 'pheatmap','dplyr','BiocManager'), repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('Seurat', 'edgeR', 'AnnotationDbi', 'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Dm.eg.db'))"

# Crea una directory di lavoro

