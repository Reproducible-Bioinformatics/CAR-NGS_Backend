#!/bin/bash

# Parametri
genome_dir="/genome/"
fastq_dir="/scratch/"
adapter_seq1="${1:-AGATCGGAAGAGCACACGTCTGAACTCCAGTCA}"
adapter_seq2="${2:-AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT}"
output_dir="/scratch/output"

# Trova il file FASTA del genoma
genome_fasta=$(find "$genome_dir" -name "*.fa" -o -name "*.fasta" | head -n 1)

if [ -z "$genome_fasta" ]; then
  echo "Error: No FASTA file found in directory $genome_dir"
  exit 1
fi

# Trova il file GTF
gtf_file=$(find "$genome_dir" -name "*.gtf" -o -name "*.gff" | head -n 1)

if [ -z "$gtf_file" ]; then
  echo "Error: No GTF/GFF file found in directory $genome_dir"
  exit 1
fi

# Crea il nome della directory di indice STAR basato sul nome del file FASTA
star_index_dir="${genome_dir}/$(basename ${genome_fasta} .fa)_star_index"

# Verifica se l'indice STAR esiste
if [ ! -d "$star_index_dir" ]; then
  echo "STAR index not found, generating index..."
  mkdir -p $star_index_dir
  STAR --runThreadN 8 \
       --runMode genomeGenerate \
       --genomeDir $star_index_dir \
       --genomeFastaFiles $genome_fasta \
       --sjdbGTFfile $gtf_file \
       --sjdbOverhang 100
  echo "STAR index generated."
else
  echo "STAR index already exists, skipping generation."
fi

# Crea la directory di output se non esiste
mkdir -p $output_dir

# Trova tutti i file FASTQ nella directory fornita
fastq_files=($(find "$fastq_dir" -name "*.fastq.gz"))

# Raggruppa i file FASTQ per sample e determina paired-end o single-end
declare -A samples
for fastq_file in "${fastq_files[@]}"; do
  # Estrai la base del nome del file (senza suffisso)
  base_name=$(basename "$fastq_file" | sed -E 's/_R?[12].*//')

  # Aggiungi il file alla lista dei campioni
  samples["$base_name"]+="$fastq_file "
done

# Inizializza la matrice di conteggi e il file di metadati
count_matrix="${output_dir}/gene_count_matrix.csv"
metadata_file="${output_dir}/Covariatesstat.csv"
echo -e "Sample,Group" > $metadata_file

header="GeneID"  # Inizializza l'intestazione della matrice

# Funzione per rilevare i gruppi basata sulle parti comuni dei nomi dei campioni
detect_group() {
  local sample=$1
  echo "$sample" | sed -E 's/[0-9_]+$//'
}

# Estrai la lunghezza del gene dal file GTF
awk '$3 == "exon" {gene_length[$10]+=$5-$4} END {for (gene in gene_length) print gene, gene_length[gene]}' $gtf_file > gene_lengths.tsv

# Processo per ogni sample
for sample_name in "${!samples[@]}"; do
  echo "Processing sample: $sample_name"

  # Crea una directory per il campione specifico
  sample_output_dir="${output_dir}/${sample_name}"
  mkdir -p $sample_output_dir

  # Ottieni i file FASTQ per il sample corrente
  fastqs=(${samples[$sample_name]})

  # Determina se sono paired-end o single-end
  if [ ${#fastqs[@]} -eq 2 ]; then
    echo "Sample $sample_name is paired-end."

    # Identifica forward e reverse
    fastq1=$(echo ${fastqs[@]} | tr ' ' '\n' | grep -E '_R?1' | head -n 1)
    fastq2=$(echo ${fastqs[@]} | tr ' ' '\n' | grep -E '_R?2' | head -n 1)
    echo "Forward read: $fastq1"
    echo "Reverse read: $fastq2"

    trimmed_fastq1="${sample_output_dir}/${sample_name}_R1_trimmed.fastq.gz"
    trimmed_fastq2="${sample_output_dir}/${sample_name}_R2_trimmed.fastq.gz"

    # Utilizza cutadapt per tagliare gli adattatori e rimuovere automaticamente le letture paired-end che diventano troppo corte
    cutadapt -j 8 -a $adapter_seq1 -A $adapter_seq2 -o $trimmed_fastq1 -p $trimmed_fastq2 --minimum-length 20 $fastq1 $fastq2

  else
    echo "Sample $sample_name is single-end."

    fastq1=${fastqs[0]}
    trimmed_fastq1="${sample_output_dir}/${sample_name}_trimmed.fastq.gz"

    cutadapt -j 8 -a $adapter_seq1 -o $trimmed_fastq1 --minimum-length 20 $fastq1
  fi

  # Allinea le letture con STAR e genera la matrice di conteggio
  if [ -n "$trimmed_fastq2" ]; then
    # Paired-end
    STAR --runThreadN 8 \
         --genomeDir $star_index_dir \
         --readFilesIn $trimmed_fastq1 $trimmed_fastq2 \
         --readFilesCommand zcat \
         --outFileNamePrefix $sample_output_dir/${sample_name}_ \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts
  else
    # Single-end
    STAR --runThreadN 8 \
         --genomeDir $star_index_dir \
         --readFilesIn $trimmed_fastq1 \
         --readFilesCommand zcat \
         --outFileNamePrefix $sample_output_dir/${sample_name}_ \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts
  fi

  # Path al file di conteggio delle espressioni geniche
  output_counts="${sample_output_dir}/${sample_name}_ReadsPerGene.out.tab"

  # Aggiungi i dati alla matrice di conteggio
  if [ "$sample_name" == "$(echo ${!samples[@]} | awk '{print $1}')" ]; then
    # Per il primo campione, aggiungi i nomi dei geni e le prime colonne
    awk 'NR>4 {print $1}' $output_counts > temp_gene_ids.tsv
    awk 'NR>4 {print $2}' $output_counts > temp_counts.tsv
    paste temp_gene_ids.tsv temp_counts.tsv | tr "\t" "," > temp_matrix.csv
    mv temp_matrix.csv $count_matrix
    header="${header},${sample_name}"  # Aggiungi il nome del campione all'intestazione
  else
    # Per gli altri campioni, aggiungi solo le colonne delle conte
    awk 'NR>4 {print $2}' $output_counts > temp_counts.tsv
    paste -d ',' $count_matrix temp_counts.tsv > temp_matrix.csv
    mv temp_matrix.csv $count_matrix
    header="${header},${sample_name}"  # Aggiungi il nome del campione all'intestazione
  fi

  rm -f temp_gene_ids.tsv temp_counts.tsv  # Rimuove i file temporanei solo se esistono

  echo "Analysis for sample $sample_name completed. Counts saved in $output_counts"

  # Rileva il gruppo
  group=$(detect_group $sample_name)

  # Scrivi i metadati nel file
  echo -e "${sample_name},${group}" >> $metadata_file

done

# Aggiungi l'intestazione alla matrice finale
sed -i "1s/^/${header}\n/" $count_matrix

echo "RNA-seq analysis pipeline completed for all samples. Count matrix saved in $count_matrix and metadata saved in $metadata_file."
