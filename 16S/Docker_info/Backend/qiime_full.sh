#!/bin/bash

# Definisci la directory dove si trovano i file fastq
FASTQ_DIR="/scratch/"

# Crea una directory per i risultati dell'allineamento
OUTPUT_DIR="${FASTQ_DIR}/aligned_results"
HTML_DIR="${OUTPUT_DIR}/html"
mkdir -p ${OUTPUT_DIR}
mkdir -p ${HTML_DIR}

# Percorso del classificatore pre-addestrato Silva 138
CLASSIFIER="/home/silva-138-99-nb-classifier.qza"

# Estensione dei file fastq da cercare
EXTENSIONS=("fastq.gz" "fastq" "fq.gz" "fq")

# Flag per verificare se sono stati trovati file
FILES_FOUND=0

# Loop attraverso ogni estensione di file
for EXT in "${EXTENSIONS[@]}"
do
    # Cerca i file R1 con l'estensione corrente
    for R1_FILE in ${FASTQ_DIR}/*_R1.${EXT}
    do
        # Verifica se esiste il file R1
        if [ ! -e "$R1_FILE" ]; then
            continue
        fi
        
        # Incrementa il flag di file trovati
        FILES_FOUND=1

        # Identifica il file R2 corrispondente
        R2_FILE="${R1_FILE/_R1/_R2}"
        
        # Controlla se il file R2 esiste
        if [ ! -e "$R2_FILE" ]; then
            echo "ATTENZIONE: File R2 corrispondente non trovato per $R1_FILE"
            continue
        fi
        
        # Estrai il nome base per i file di output
        BASE_NAME=$(basename ${R1_FILE} _R1.${EXT})
        
        # Crea una directory temporanea per il manifest
        TEMP_DIR="${FASTQ_DIR}/temp_${BASE_NAME}"
        mkdir -p ${TEMP_DIR}
        
        # Crea il manifest file (TSV)
        MANIFEST_FILE="${TEMP_DIR}/manifest.tsv"
        echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > ${MANIFEST_FILE}
        echo -e "${BASE_NAME}\t${R1_FILE}\t${R2_FILE}" >> ${MANIFEST_FILE}
        
        # Esegui l'importazione dei dati
        IMPORT_CMD="qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ${MANIFEST_FILE} --output-path ${OUTPUT_DIR}/${BASE_NAME}-paired-end-demux.qza --input-format PairedEndFastqManifestPhred33V2"
        echo "Eseguendo: ${IMPORT_CMD}"
        eval ${IMPORT_CMD}

        # Riassunto del demux
        DEMUX_CMD="qiime demux summarize --i-data ${OUTPUT_DIR}/${BASE_NAME}-paired-end-demux.qza --o-visualization ${OUTPUT_DIR}/${BASE_NAME}-paired-end-demux.qzv"
        echo "Eseguendo: ${DEMUX_CMD}"
        eval ${DEMUX_CMD}

        # Denoising con DADA2
        DENOISE_CMD="qiime dada2 denoise-paired --i-demultiplexed-seqs ${OUTPUT_DIR}/${BASE_NAME}-paired-end-demux.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 250 --p-trunc-len-r 250 --o-table ${OUTPUT_DIR}/${BASE_NAME}-table.qza --o-representative-sequences ${OUTPUT_DIR}/${BASE_NAME}-rep-seqs.qza --o-denoising-stats ${OUTPUT_DIR}/${BASE_NAME}-denoising-stats.qza"
        echo "Eseguendo: ${DENOISE_CMD}"
        eval ${DENOISE_CMD}
        
        # Tabulazione dei dati di denoising
        TABULATE_CMD="qiime metadata tabulate --m-input-file ${OUTPUT_DIR}/${BASE_NAME}-denoising-stats.qza --o-visualization ${OUTPUT_DIR}/${BASE_NAME}-denoising-stats.qzv"
        echo "Eseguendo: ${TABULATE_CMD}"
        eval ${TABULATE_CMD}

        # Classificazione tassonomica
        TAXONOMY_CMD="qiime feature-classifier classify-sklearn --i-classifier ${CLASSIFIER} --i-reads ${OUTPUT_DIR}/${BASE_NAME}-rep-seqs.qza --o-classification ${OUTPUT_DIR}/${BASE_NAME}-taxonomy.qza"
        echo "Eseguendo: ${TAXONOMY_CMD}"
        eval ${TAXONOMY_CMD}

        # Visualizza i risultati della tassonomia (senza file di metadata)
        TAXA_BARPLOT_CMD="qiime taxa barplot --i-table ${OUTPUT_DIR}/${BASE_NAME}-table.qza --i-taxonomy ${OUTPUT_DIR}/${BASE_NAME}-taxonomy.qza --o-visualization ${OUTPUT_DIR}/${BASE_NAME}-taxa-bar-plots.qzv"
        echo "Eseguendo: ${TAXA_BARPLOT_CMD}"
        eval ${TAXA_BARPLOT_CMD}

        # Esporta il contenuto del file qzv in una directory HTML
        EXPORT_CMD="qiime tools export --input-path ${OUTPUT_DIR}/${BASE_NAME}-taxa-bar-plots.qzv --output-path ${HTML_DIR}/${BASE_NAME}-taxa-bar-plots"
        echo "Eseguendo: ${EXPORT_CMD}"
        eval ${EXPORT_CMD}

        # Pulisce la directory temporanea
        rm -r ${TEMP_DIR}

    done
done

# Controlla se non sono stati trovati file
if [ $FILES_FOUND -eq 0 ]; then
    echo "ATTENZIONE: Nessun file FASTQ trovato nella directory ${FASTQ_DIR} con le estensioni ${EXTENSIONS[*]}"
fi
