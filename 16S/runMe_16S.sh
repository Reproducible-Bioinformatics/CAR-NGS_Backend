#!/bin/bash

# Definisci le directory locali
LOCAL_FASTQ_DIR=$1   # Directory dove ci sono i dati. I dati DEVONO essere pair end fastq.gz, fastq, fq.gz o fq. Il nome dei files deve avere la seguente forma : nome_R1.fastq.gz oppure nome_R2.fastq.gz
# Esegui il container Docker con le directory montate
docker run --rm -it \
-v ${LOCAL_FASTQ_DIR}:/scratch \
--user $(id -u):$(id -g) \
repbioinfo/qiime2023 /home/qiime_full.sh
