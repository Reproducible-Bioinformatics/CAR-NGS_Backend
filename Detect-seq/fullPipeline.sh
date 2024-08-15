#!/bin/bash

# Trova il numero massimo degli esperimenti esistenti e incrementa di 1
exp_num=1
for dir in "$(pwd)"/Exp*; do
  if [[ -d "$dir" && "$dir" =~ Exp([0-9]+) ]]; then
    num=${BASH_REMATCH[1]}
    if (( num > exp_num )); then
      exp_num=$((num + 1))
    fi
  fi
done

# Crea la directory necessaria per l'esperimento
exp_dir="$(pwd)/Exp${exp_num}"
mkdir -p "${exp_dir}"

# Definisce le variabili per il threshold e gli adattatori
genome_dir=$1
threshold=$3
adapt1=$4
adapt2=$5

# Esegui il container Docker, montando le cartelle corrette
docker run --rm -i \
  -v "${genome_dir}:/genome" \
  -v "${exp_dir}":/scratch \
  -v "$2":/scratch/raw.fastq:ro \
  repbioinfo/detectseq  \
  /home/detectSeq.sh $threshold $adapt1 $adapt2
