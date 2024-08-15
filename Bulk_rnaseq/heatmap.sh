#!/bin/bash

# Funzione per generare la heatmap con Docker
esegui_heatmap_docker() {
    # Determina il percorso della directory corrente
    CURRENT_DIR=$(pwd)

    # Directory dei dati
    DATA_DIR=$(realpath "$CURRENT_DIR/Data")

    # Controlla se la directory dei dati esiste
    if [ ! -d "$DATA_DIR" ]; then
      echo "Error: Data directory ($DATA_DIR) not found."
      exit 1
    fi

    # Nome dell'immagine Docker
    DOCKER_IMAGE="repbioinfo/rnaseqstar_v2"

    # Comando per eseguire Docker con le cartelle montate correttamente
    docker run --rm \
        -v "$DATA_DIR:/scratch" \
        "$DOCKER_IMAGE" \
        Rscript /home/heatmap.R gene_count_matrix.csv Covariatesstat.csv
}

# Esempio di utilizzo della funzione
esegui_heatmap_docker

