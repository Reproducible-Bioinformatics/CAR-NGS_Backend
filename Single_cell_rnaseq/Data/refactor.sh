#!/bin/bash

# Cartelle di input e output
fastq_dir="/scratch/"
fixed_fastq_dir="/scratch/fix.fastq"
log_file="/scratch/fix.fastq/rename_log.txt"

# Creazione della cartella di output se non esiste
mkdir -p "$fixed_fastq_dir"

# Inizializzazione del file di log
echo "Log delle trasformazioni FASTQ" > "$log_file"
echo "Data: $(date)" >> "$log_file"
echo "-----------------------------" >> "$log_file"

# Funzione per estrarre l'estensione del file e normalizzare
normalize_extension() {
  if [[ "$1" == *.fq.gz ]]; then
    echo ".fastq.gz"
  elif [[ "$1" == *.fastq.gz ]]; then
    echo ".fastq.gz"
  fi
}

# Funzione per controllare e correggere i nomi dei file
fix_fastq_name() {
  fastq_file="$1"
  ext=$(normalize_extension "$fastq_file")
  base_name=$(basename "$fastq_file" | sed "s/.fq.gz$//" | sed "s/.fastq.gz$//")

  # Estrai il nome del campione
  sample_name=$(echo "$base_name" | awk -F'_' '{print $1}')

  # Controlla se esiste "_S" seguito da un numero
  if [[ "$base_name" =~ _S[0-9]+ ]]; then
    s_index=$(echo "$base_name" | grep -o '_S[0-9]\+')
  else
    s_index="_S1"
  fi

  # Controlla se esiste "_L" seguito da un numero
  if [[ "$base_name" =~ _L[0-9]{3} ]]; then
    lane=$(echo "$base_name" | grep -o '_L[0-9]\{3\}')
  else
    lane="_L001"
  fi

  # Controlla se esiste "_R" seguito da un numero, oppure "_1" o "_2"
  if [[ "$base_name" =~ _R[12] ]]; then
    read_type=$(echo "$base_name" | grep -o '_R[12]')
  elif [[ "$base_name" =~ _1$ ]]; then
    read_type="_R1"
  elif [[ "$base_name" =~ _2$ ]]; then
    read_type="_R2"
  else
    # Se non si trova né "_R1/_R2" né "_1/_2", usa "_R1" come default
    read_type="_R1"
  fi

  # Controlla se esiste "_001" alla fine
  if [[ "$base_name" =~ _[0-9]{3}$ ]]; then
    chunk=$(echo "$base_name" | grep -o '_[0-9]\{3\}$')
  else
    chunk="_001"
  fi

  # Costruisci il nuovo nome del file con l'estensione normalizzata
  new_base_name="${sample_name}${s_index}${lane}${read_type}${chunk}.fastq.gz"
  new_fastq_file="${fixed_fastq_dir}/${new_base_name}"

  # Verifica se il link simbolico esiste già
  if [ -e "$new_fastq_file" ]; then
    echo "Warning: $new_fastq_file already exists, skipping..."
  else
    # Crea un link simbolico al file con il nuovo nome
    ln -s "$fastq_file" "$new_fastq_file"
    echo "Created symlink from $fastq_file to $new_fastq_file"

    # Registra la trasformazione nel file di log
    echo "File: $fastq_file -> $new_fastq_file" >> "$log_file"
  fi
}

# Trova tutti i file FASTQ (fastq.gz e fq.gz) e correggi i nomi
fastq_files=($(find "$fastq_dir" -name "*.fastq.gz" -o -name "*.fq.gz"))

for fastq_file in "${fastq_files[@]}"; do
  fix_fastq_name "$fastq_file"
done

echo "Tutti i file FASTQ sono stati controllati e corretti, se necessario. I link simbolici sono stati creati."
echo "Il log delle trasformazioni è stato salvato in $log_file."
