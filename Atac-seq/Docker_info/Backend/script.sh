#!/bin/bash

# Parametri
DATA_DIR="/scratch" # Directory contenente i file FASTQ
GENOME_DIR="/genomes" # Directory contenente il file FASTA del genoma
OUTPUT_DIR="/scratch/results" # Directory di output fissa
THREADS=${1:-8} # Numero di thread da utilizzare

# Trova il file FASTA del genoma nella directory /genomes
GENOME_FASTA=$(find $GENOME_DIR -name "*.fa" -o -name "*.fasta" | head -n 1)
if [ -z "$GENOME_FASTA" ]; then
    echo "Errore: File FASTA del genoma non trovato nella directory /genomes."
    exit 1
fi

# Creare directory per output
mkdir -p "$OUTPUT_DIR"
QUALITY_DIR="$OUTPUT_DIR/Quality_ATAC"
mkdir -p "$QUALITY_DIR"
PEAK_OUTPUT_DIR="$OUTPUT_DIR/peaks"
mkdir -p "$PEAK_OUTPUT_DIR"

# Step 1: Controllo della presenza dell'indice del genoma
GENOME_INDEX_DIR="$GENOME_DIR/index"
if [ ! -f "$GENOME_INDEX_DIR/genome_index.1.bt2" ]; then
    echo "Indice del genoma non trovato. Creazione in corso..."
    mkdir -p "$GENOME_INDEX_DIR"
    bowtie2-build "$GENOME_FASTA" "$GENOME_INDEX_DIR/genome_index"
    echo "Indice del genoma creato."
else
    echo "Indice del genoma trovato. Procedo con l'analisi."
fi

# Step 2: Rileva file FASTQ nella cartella /scratch
FASTQ_FILES=($(ls "$DATA_DIR"/*.fastq.gz))
PAIR_COUNT=${#FASTQ_FILES[@]}

if [ $PAIR_COUNT -eq 0 ]; then
    echo "Errore: Nessun file FASTQ trovato."
    exit 1
fi

# Step 3: Controllo di qualità
for FILE in "${FASTQ_FILES[@]}"; do
    fastqc -o "$QUALITY_DIR" "$FILE"
done

# Step 4: Allineamento e generazione BAM file con soft clipping
for ((i=0; i<PAIR_COUNT; i++)); do
    FASTQ_FILE="${FASTQ_FILES[i]}"
    SAMPLE_NAME=$(basename "$FASTQ_FILE" | cut -d'_' -f1)
    BAM_FILE="$OUTPUT_DIR/${SAMPLE_NAME}.sorted.bam"
    NOORG_BAM_FILE="$OUTPUT_DIR/${SAMPLE_NAME}.sorted.noorg.bam"

    # Verifica se esiste un file paired-end
    PAIR_FILE=$(echo "$FASTQ_FILE" | sed 's/_R1_/_R2_/')
    if [ -f "$PAIR_FILE" ]; then
        # Bowtie2 Alignment paired-end
        bowtie2 --local --threads "$THREADS" -x "$GENOME_INDEX_DIR/genome_index" \
            -1 "$FASTQ_FILE" -2 "$PAIR_FILE" | \
        samtools view --threads "$THREADS" -bS - | \
        samtools sort --threads "$THREADS" -o "$BAM_FILE"
    else
        # Bowtie2 Alignment single-end
        bowtie2 --local --threads "$THREADS" -x "$GENOME_INDEX_DIR/genome_index" \
            -U "$FASTQ_FILE" | \
        samtools view --threads "$THREADS" -bS - | \
        samtools sort --threads "$THREADS" -o "$BAM_FILE"
    fi

    samtools index "$BAM_FILE"

    # Rimuovi letture mitocondriali e cloroplastidiche
    samtools idxstats "$BAM_FILE" | cut -f1 | grep -v Mt | grep -v Pt | xargs samtools view --threads "$THREADS" -b "$BAM_FILE" > "$NOORG_BAM_FILE"
    samtools index "$NOORG_BAM_FILE"
done

# Step 5: Chiamata dei picchi
for ((i=0; i<PAIR_COUNT; i++)); do
    FASTQ_FILE="${FASTQ_FILES[i]}"
    SAMPLE_NAME=$(basename "$FASTQ_FILE" | cut -d'_' -f1)
    NOORG_BAM_FILE="$OUTPUT_DIR/${SAMPLE_NAME}.sorted.noorg.bam"

    PAIR_FILE=$(echo "$FASTQ_FILE" | sed 's/_R1_/_R2_/')
    if [ -f "$PAIR_FILE" ]; then
        # Chiamata dei picchi paired-end
        macs2 callpeak -t "$NOORG_BAM_FILE" -q 0.05 --broad -f BAMPE -n "$SAMPLE_NAME" -B --trackline --outdir "$PEAK_OUTPUT_DIR"
    else
        # Chiamata dei picchi single-end
        macs2 callpeak -t "$NOORG_BAM_FILE" -q 0.05 --broad -f BAM -n "$SAMPLE_NAME" -B --trackline --outdir "$PEAK_OUTPUT_DIR"
    fi
done

# Step 6: Conversione dei file bedgraph in bigwig (opzionale)
for ((i=0; i<PAIR_COUNT; i++)); do
    SAMPLE_NAME=$(basename "${FASTQ_FILES[i]}" | cut -d'_' -f1)
    BEDGRAPH_FILE="$PEAK_OUTPUT_DIR/${SAMPLE_NAME}_treat_pileup.bdg"
    CLIPPED_BEDGRAPH_FILE="$PEAK_OUTPUT_DIR/${SAMPLE_NAME}_treat_pileup.clipped.sorted.bdg"
    BIGWIG_FILE="$PEAK_OUTPUT_DIR/${SAMPLE_NAME}_treat_pileup.clipped.sorted.bw"
    CHR_SIZES_FILE="$OUTPUT_DIR/At_chr.sizes"

    bioawk -c fastx '{print $name, length($seq)}' "$GENOME_FASTA" > "$CHR_SIZES_FILE"
    bedtools slop -i "$BEDGRAPH_FILE" -g "$CHR_SIZES_FILE" -b 0 | bedClip stdin "$CHR_SIZES_FILE" "$CLIPPED_BEDGRAPH_FILE"
    sort -k1,1 -k2,2n "$CLIPPED_BEDGRAPH_FILE" > "$CLIPPED_BEDGRAPH_FILE.sorted"
    bedGraphToBigWig "$CLIPPED_BEDGRAPH_FILE.sorted" "$CHR_SIZES_FILE" "$BIGWIG_FILE"
done

echo "Analisi ATAC-seq completata!"
