# Usa un'immagine base di Ubuntu
FROM ubuntu:22.04

# Mantieni l'immagine aggiornata
RUN apt-get update && apt-get upgrade -y

# Installa le dipendenze necessarie
RUN apt-get install -y \
    wget \
    curl \
    unzip \
    git \
    build-essential \
    python3 \
    python3-pip \
    openjdk-11-jre \
    samtools \
    bowtie2 \
    bedtools \
    bioawk \
    zlib1g-dev

# Installa FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    mv FastQC /opt/ && \
    chmod +x /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc

# Installa MACS2 (richiede Python e pip)
RUN pip3 install macs2

# Installa UCSC bedClip e bedGraphToBigWig
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedClip && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && \
    chmod +x bedClip bedGraphToBigWig && \
    mv bedClip /usr/local/bin/ && \
    mv bedGraphToBigWig /usr/local/bin/

# Imposta la directory di lavoro
WORKDIR /data

# Comando di default quando si avvia il container
CMD ["bash"]
