#!/bin/bash

# Ricombina tutti i file decompressi in un unico archivio
cat Mus_musculus.GRCm39.dna_rm.toplevel.tar.gz.part* > Mus_musculus.GRCm39.dna_rm.toplevel.tar.gz

# Rimuovi le parti decompressi per pulizia
rm Mus_musculus.GRCm39.dna_rm.toplevel.tar.gz.part*

# Estrai l'archivio tar.gz ricombinato
tar -xzvf Mus_musculus.GRCm39.dna_rm.toplevel.tar.gz

# Rimuovi l'archivio ricombinato
rm Mus_musculus.GRCm39.dna_rm.toplevel.tar.gz

echo "File ricombinato e estratto in Mus_musculus.GRCm39.dna_rm.toplevel.fa"
