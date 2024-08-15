#!/bin/bash

# Funzione per eseguire l'allineamento e l'indicizzazione con Docker
esegui_index_align_docker() {
    # Determina il percorso della directory corrente
    CURRENT_DIR=$(pwd)

    # Directory dei dati e del genoma
    DATA_DIR=$(realpath "$CURRENT_DIR/Data")
    GENOME_DIR=$(realpath "$CURRENT_DIR/Genome")

    # Controlla se le directory esistono
    if [ ! -d "$DATA_DIR" ]; then
      echo "Error: Data directory ($DATA_DIR) not found."
      exit 1
    fi

    if [ ! -d "$GENOME_DIR" ]; then
      echo "Error: Genome directory ($GENOME_DIR) not found."
      exit 1
    fi

    # Nome dell'immagine Docker
    DOCKER_IMAGE="repbioinfo/rnaseqstar_v2"

    # Comando per eseguire Docker con le cartelle montate correttamente
    docker run --rm \
        -v "$DATA_DIR:/scratch" \
        -v "$GENOME_DIR:/genome" \
        "$DOCKER_IMAGE" \
        bash /home/index_align.sh
}

# Esempio di utilizzo della funzione
esegui_index_align_docker
