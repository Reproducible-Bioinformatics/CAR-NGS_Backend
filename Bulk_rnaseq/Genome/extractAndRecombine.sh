 
gunzip part_*.gz
cat part_* > Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa
#!/bin/bash
rm part_*

# Decomprimi tutti i file con il prefisso specificato
gunzip split_gtf_*.gz

# Ricombina tutti i file decompressi in un unico file
cat split_gtf_* > Drosophila_melanogaster.BDGP6.46.112.gtf

# Rimuovi le parti decompressi per pulizia
rm split_gtf_*

echo "File ricombinato in Drosophila_melanogaster.BDGP6.46.112.gtf"
