#!/bin/bash

# Funzione per eseguire l'analisi ATAC-seq tramite Docker
esegui_atac_seq_docker() {
    # Parametri di input
    FASTQ_DIR=$1    # Directory contenente i file FASTQ
    FASTA_DIR=$2    # Directory contenente il file FASTA del genoma
    THREADS=${3:-8} # Numero di thread da utilizzare (opzionale, default: 8)

    # Nome dell'immagine Docker
    DOCKER_IMAGE="repbioinfo/atacseq" # Nome dell'immagine Docker

    # Comando per eseguire Docker con le cartelle montate correttamente
    docker run --rm \
        -v "$FASTQ_DIR:/scratch" \
        -v "$FASTA_DIR:/genomes" \
        -v "$FASTQ_DIR/results:/scratch/results" \
        "$DOCKER_IMAGE" \
        bash -c "/home/script.sh $THREADS"
}

# Esempio di utilizzo della funzione
# esegui_atac_seq_docker /path/to/fastq_files /path/to/fasta_files 8
