#!/bin/bash

WORKING_DIR=/data/scratch/
cd $WORKING_DIR

SCRIPTS_DIR=/home/jspacker/sci-RNA-seq-pipeline-scripts
DATAMASH_PATH=/home/jspacker/datamash/datamash

STAR_INDEX=/home/GENOMES
GENE_MODEL_DIR=/home/GENOMES

RT_BARCODES_FILE=$WORKING_DIR/sci-RNA-seq-8.RT.oligos
SAMPLE_NAME="$1"
echo "$SAMPLE_NAME" >$WORKING_DIR/combinatorial.indexing.key
BATCH_SIZE=1
#ILLUMINA_DIR=/home/220125_NB501050_0349_AHLWKVBGXK
UMI_PER_CELL_CUTOFF=$2

# #-------------------------------------------------------------------------------
# # Extract FASTQ files from BCL files
# #calculate reads loss
# #-------------------------------------------------------------------------------
# sudo apt-get install alien --assume-yes
# sudo alien bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm
# sudo dpkg -i bcl2fastq2_0v2.20.0.422-2_amd64.deb
#
# cd $WORKING_DIR
# rm -r fastq
# mkdir fastq
#
# sh $SCRIPTS_DIR/run-bcl2fastq.sh \
#     $ILLUMINA_DIR/ \
#     $WORKING_DIR/fastq \
#     $WORKING_DIR/SampleSheet.csv

#-------------------------------------------------------------------------------
# Put read 1 info (RT well, UMI) into read 2 read name
# NB! fastq file name must start with NAME OF THE WELL e.g. A01 (see script at the end to modify the names)
#-------------------------------------------------------------------------------
#yes Y | apt-get install gawk

cd $WORKING_DIR
#rm -r combined-fastq
#rm -r file-lists-for-r1-info-munging
#rm -r put-r1-info-in-r2-logs
mkdir combined-fastq
mkdir file-lists-for-r1-info-munging
mkdir put-r1-info-in-r2-logs

ls fastq/ | grep _R1_ | grep -v Undetermined | split -l $BATCH_SIZE -d - file-lists-for-r1-info-munging/

ls file-lists-for-r1-info-munging | while read BATCH
do
    bash $SCRIPTS_DIR/put-read1-info-in-read2.sh $WORKING_DIR/fastq $WORKING_DIR/file-lists-for-r1-info-munging/$BATCH $SCRIPTS_DIR/ $RT_BARCODES_FILE $WORKING_DIR/combinatorial.indexing.key $WORKING_DIR/combined-fastq $WORKING_DIR/put-r1-info-in-r2-logs

done

#rm ./RT_barcode_stats.txt
ls file-lists-for-r1-info-munging | while read BATCH
do
    paste file-lists-for-r1-info-munging/$BATCH put-r1-info-in-r2-logs/$BATCH >> ./RT_barcode_stats.txt
done

#-------------------------------------------------------------------------------
# Trim poly-A tails
#docker pull quay.io/biocontainers/trim-galore:0.4.1--pl5.22.0_0
#-------------------------------------------------------------------------------
#curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz
#tar xvzf trim_galore.tar.gz --no-same-owner
#yes Y | sudo apt install cutadapt


#rm -r file-lists-for-trimming
#rm -r trimmed-fastq
mkdir file-lists-for-trimming
mkdir trimmed-fastq

ls combined-fastq/ | split -l $BATCH_SIZE -d - file-lists-for-trimming/

ls file-lists-for-trimming | while read BATCH; do
    sh $SCRIPTS_DIR/run-trim-galore.sh     \
        $WORKING_DIR/combined-fastq                         \
        $WORKING_DIR/file-lists-for-trimming/$BATCH         \
        $WORKING_DIR/trimmed-fastq
done

#rm ./trimming_report.txt
tail -n 3 ./trimmed-fastq/*report.txt >> trimming_report.txt

#-------------------------------------------------------------------------------
# Align reads using STAR
#-------------------------------------------------------------------------------
# wget https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz
# tar -xzf 2.7.10a.tar.gz --no-same-owner
# cd STAR-2.7.10a
# cd ./source
# make STAR
# 
# cd $WORKING_DIR

#install samtools
# wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2
# tar xvjf samtools-1.3.tar.bz2 --no-same-owner
# cd samtools-1.3
# ./configure
# make
# sudo make install
export PATH="$PATH:/home/samtools-1.3"
# cd $WORKING_DIR
#generate gtf and fasta for human and mouse
#cd GENOMES
#awk '/^>/{$1=$1"M"} 1' Mus_musculus.GRCm39.dna.primary_assembly.fa >> mouse_no_under.fa
#awk '{$1=$1"M"} 1' Mus_musculus.GRCm39.109.gtf >> mouse_no_under_gtf.gtf

#generate star index from gtf and fasta, only if you want to use a new STAR index
#sudo apt-get install -y rsem
#THREADS=12
#RSEM="/home/RSEM-1.2.30"
#chmod 777 /home/RSEM-1.2.30/*
#cd /home/star
#$RSEM/rsem-prepare-reference -p $THREADS --star --star-path /home/STAR-2.7.10a/source --gtf /home/star/Homo_sapiens.GRCh38.109.gtf /home/star/Homo_sapiens.GRCh38.dna.primary_assembly.fa
#/home/STAR-2.7.10a/source/STAR --runMode genomeGenerate --runThreadN 12 --genomeDir /home/star  --genomeFastaFiles /home/star/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --outTmpDir /home/star_new --sjdbGTFfile /home/star/Homo_sapiens.GRCh38.109.gtf
#$RSEM/rsem-prepare-reference -p $THREADS --star --star-path /home/STAR-2.5.1b/source --gtf /home/GENOMES/Homo_sapiens.GRCh38.109.gtf /home/GENOMES/Homo_sapiens.GRCh38.dna.primary_assembly.fa /home/star

#alignin with STAR
# wget https://github.com/alexdobin/STAR/archive/2.5.1b.tar.gz
# tar -xzf 2.5.1b.tar.gz --no-same-owner
# cd STAR-2.5.1b
# cd ./source
# make STAR

cd $WORKING_DIR

#rm -r aligned-reads
mkdir aligned-reads

sh $SCRIPTS_DIR/STAR-alignReads-moreMem.sh \
    $WORKING_DIR/trimmed-fastq                              \
    $STAR_INDEX                                           \
    $WORKING_DIR/aligned-reads

# comprehensive statistics on alignment
#rm ./alignment_report.txt
#rm ./alignment_report_complete.txt
#rm ./sample_list.txt
cat ./aligned-reads/*final.out | grep "Uniquely mapped reads %" >> ./alignment_report.txt
ls file-lists-for-r1-info-munging | while read BATCH
do
    cat file-lists-for-r1-info-munging/$BATCH >> ./sample_list.txt
done

paste ./sample_list.txt ./alignment_report.txt >> ./alignment_report_complete.txt

#-------------------------------------------------------------------------------
# Filter ambiguously-mapped reads and sort BAM files
# Also count rRNA reads
#-------------------------------------------------------------------------------

cd $WORKING_DIR

#rm -r file-lists-for-samtools-sort
#rm -r aligned-reads-filtered-sorted
#rm -r rRNA-read-counts
mkdir file-lists-for-samtools-sort
mkdir aligned-reads-filtered-sorted
mkdir rRNA-read-counts

#install samtools
# wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2
# tar xvjf samtools-1.7.tar.bz2 --no-same-owner
# cd samtools-1.7
# ./configure
# make
# sudo make install
export PATH="$PATH:/data/scratch/samtools-1.7"
# cd $WORKING_DIR

ls aligned-reads/ | grep "[.]Aligned[.]out[.]bam$" | split -l $BATCH_SIZE -d - file-lists-for-samtools-sort/

ls file-lists-for-samtools-sort | while read BATCH; do
    $SCRIPTS_DIR/samtools-filter-sort.sh    \
        $WORKING_DIR/aligned-reads                              \
        $WORKING_DIR/file-lists-for-samtools-sort/$BATCH        \
        $WORKING_DIR/aligned-reads-filtered-sorted
done

#install bedtools
# sudo apt install python-is-python3
# wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz
# tar -zxvf bedtools-2.27.1.tar.gz --no-same-owner
# cd bedtools2
# make

#install bedops to generate bed
#sudo apt-get install -y bedops

#generate bed
#gtf2bed < $GENE_MODEL_DIR/Homo_sapiens.GRCh38.101.gtf > v43.bed
#awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' $GENE_MODEL_DIR/Homo_sapiens.GRCh38.101.gtf | gtf2bed - > $GENE_MODEL_DIR/output.bed

cd $WORKING_DIR

ls file-lists-for-samtools-sort | while read BATCH; do
    $SCRIPTS_DIR/count-rRNA-reads.sh        \
        $WORKING_DIR/aligned-reads                              \
        $WORKING_DIR/file-lists-for-samtools-sort/$BATCH        \
        $GENE_MODEL_DIR/latest.rRNA.gene.regions.union.bed     # \
        #"$WORKING_DIR/rRNA-read-counts"
done

#stats rRNA
#rm ./rRNA_report.txt
#rm ./rRNA_report_complete.txt
cat $WORKING_DIR/rRNA-read-counts/* >> ./rRNA_report.txt
paste ./sample_list.txt ./rRNA_report.txt >> ./rRNA_report_complete.txt

cat ./rRNA_report.txt
#-------------------------------------------------------------------------------
# Split reads in BAM files into BED intervals 6
#-------------------------------------------------------------------------------

cd $WORKING_DIR

#rm -r file-lists-for-rmdup
#rm -r aligned-reads-rmdup-split-bed
mkdir file-lists-for-rmdup
mkdir aligned-reads-rmdup-split-bed

ls aligned-reads-filtered-sorted/ | grep "[.]bam$" | split -l $BATCH_SIZE -d - file-lists-for-rmdup/

ls file-lists-for-rmdup | while read BATCH; do
    $SCRIPTS_DIR/rmdup-and-make-split-bed.sh    \
        $WORKING_DIR/aligned-reads-filtered-sorted                  \
        $WORKING_DIR/file-lists-for-rmdup/$BATCH                    \
        $SCRIPTS_DIR/                                               \
        $WORKING_DIR/aligned-reads-rmdup-split-bed
done


#-------------------------------------------------------------------------------
# Assign reads to genes, using the BED files as input 7
#-------------------------------------------------------------------------------

cd $WORKING_DIR

#rm -r file-lists-for-assign-reads-to-genes
#rm -r unique-read-to-gene-assignments
mkdir file-lists-for-assign-reads-to-genes
mkdir unique-read-to-gene-assignments

ls aligned-reads-rmdup-split-bed/ | grep "[.]bed$" | split -l $BATCH_SIZE -d - file-lists-for-assign-reads-to-genes/

ls file-lists-for-assign-reads-to-genes | while read BATCH; do
    $SCRIPTS_DIR/assign-reads-to-genes.sh       \
        $WORKING_DIR/aligned-reads-rmdup-split-bed                  \
        $WORKING_DIR/file-lists-for-assign-reads-to-genes/$BATCH    \
        $GENE_MODEL_DIR/latest.exons.bed                            \
        $GENE_MODEL_DIR/latest.genes.bed                            \
        $SCRIPTS_DIR/                                               \
        $WORKING_DIR/unique-read-to-gene-assignments
done


#-------------------------------------------------------------------------------
# Compute the duplication rate and proportion of reads that are from rRNA 8
#-------------------------------------------------------------------------------

cd $WORKING_DIR

#rm -r UMI-counts-by-sample
#rm -r file-lists-for-UMI-counting
mkdir UMI-counts-by-sample
mkdir file-lists-for-UMI-counting

ls aligned-reads-rmdup-split-bed/ | while read FILE; do
    PCR_WELL=`basename $FILE .bed`
    echo "$PCR_WELL"
done \
| split -l $BATCH_SIZE -d - file-lists-for-UMI-counting/

ls file-lists-for-UMI-counting | while read BATCH; do
    $SCRIPTS_DIR/count-UMI-per-sample.sh    \
        $WORKING_DIR/aligned-reads-filtered-sorted              \
        $WORKING_DIR/aligned-reads-rmdup-split-bed              \
        $WORKING_DIR/file-lists-for-UMI-counting/$BATCH         \
        $WORKING_DIR/UMI-counts-by-sample
done

#stats UMI
#rm ./UMI_report.txt
#rm ./UMI_report_complete.txt
cat $WORKING_DIR/UMI-counts-by-sample/*UMI.count >> ./UMI_report.txt
paste ./sample_list.txt ./UMI_report.txt >> ./UMI_report_complete.txt

#stats reads
#rm ./reads_report_complete.txt
for file in ./UMI-counts-by-sample/*; do
    # Check if it's a regular file
    if [ -f "$file" ]; then
        # Print the file name to the output file
        echo "File: $file" >> ./reads_report_complete.txt

        # Print the content of the file to the output file
        cat "$file" >> ./reads_report_complete.txt
    fi
done

#-------------------------------------------------------------------------------

cat UMI-counts-by-sample/*.UMI.count | sort -k1,1 \
| $DATAMASH_PATH -g 1 sum 2 \
>total.UMI.count.by.sample

cat UMI-counts-by-sample/*.read.count | sort -k1,1 \
| $DATAMASH_PATH -g 1 sum 2 \
>total.read.count.by.sample

#rm -r final-output
mkdir final-output

cat rRNA-read-counts/* | sort -k1,1 | $DATAMASH_PATH -g 1 sum 2 sum 3 \
| join - total.UMI.count.by.sample \
| join - total.read.count.by.sample \
| awk 'BEGIN {
    printf "%-18s    %11s    %8s    %10s    %8s\n",
        "sample", "n.reads", "pct.rRNA", "n.UMI", "dup.rate";
} {
    printf "%-18s    %11d    %7.1f%%    %10d    %7.1f%%\n",
        $1, $3, 100 * $2/$3, $4, 100 * (1 - $4/$5);
}' \
>final-output/rRNA.and.dup.rate.stats

cat final-output/rRNA.and.dup.rate.stats

#-------------------------------------------------------------------------------
# Make knee plots 9
#-------------------------------------------------------------------------------

cd $WORKING_DIR

#
# will take ~5 minutes per ~25M UMIs reported in final-output/rRNA.and.dup.rate.stats
#
#yes Y | apt-get install autoconf
#sudo apt-get install gawk

cat unique-read-to-gene-assignments/* | awk '$3 == "exonic" || $3 == "intronic" {
    split($1, arr, "|");
    printf "%s|%s_%s_%s\n", arr[2], arr[3], arr[4], arr[5];
}' \
| sort -k1,1 -S 8G --compress-program=/bin/gzip \
| $DATAMASH_PATH -g 1 count 1 \
| tr '|' '\t' \
>UMIs.per.cell.barcode

#rm -r final-output/knee-plots
mkdir final-output/knee-plots

# sudo echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" | sudo tee -a /etc/apt/sources.list
# gpg --keyserver keyserver.ubuntu.com --recv-key E298A3A825C0D65DFD57CBB651716619E084DAB9
# gpg -a --export E298A3A825C0D65DFD57CBB651716619E084DAB9 | sudo apt-key add -
# sudo apt-get update
# sudo apt-get install -y r-base r-base-dev r-cran-xml r-cran-rjava libcurl4-openssl-dev
# sudo apt-get install -y libssl-dev libxml2-dev openjdk-7-* libgdal-dev libproj-dev libgsl-dev
# sudo apt-get install -y xml2 default-jre default-jdk mesa-common-dev libglu1-mesa-dev freeglut3-dev
# sudo apt-get install -y mesa-common-dev libx11-dev r-cran-rgl r-cran-rglpk r-cran-rsymphony r-cran-plyr
# sudo apt-get install -y  r-cran-reshape  r-cran-reshape2 r-cran-rmysql
# sudo R CMD javareconf
# Rscript -e 'install.packages("dplyr", repos="https://cloud.r-project.org")'
# Rscript -e 'install.packages("ggplot2", repos="https://cloud.r-project.org")'

Rscript $SCRIPTS_DIR/knee-plot.R            \
    UMIs.per.cell.barcode                   \
    $WORKING_DIR/final-output/knee-plots

cd final-output && tar czf ../knee-plots.tar.gz knee-plots/ && cd ..

#UMI_PER_CELL_CUTOFF=500

Rscript $SCRIPTS_DIR/knee-plot.R            \
    UMIs.per.cell.barcode                   \
    $WORKING_DIR/final-output/knee-plots    \
    $UMI_PER_CELL_CUTOFF

cd final-output && tar czf ../knee-plots.tar.gz knee-plots/ && cd ..


#-------------------------------------------------------------------------------
# Make the final UMI count matrix
#-------------------------------------------------------------------------------

cp $GENE_MODEL_DIR/gene.annotations final-output/gene.annotations

#rm -r file-lists-for-UMI-count-rollup
#rm -r UMI-count-rollup
mkdir file-lists-for-UMI-count-rollup
mkdir UMI-count-rollup

ls unique-read-to-gene-assignments/ | split -l $BATCH_SIZE -d - file-lists-for-UMI-count-rollup/

ls file-lists-for-UMI-count-rollup | while read BATCH; do
    $SCRIPTS_DIR/UMI-count-rollup.sh        \
        $WORKING_DIR/unique-read-to-gene-assignments            \
        $WORKING_DIR/file-lists-for-UMI-count-rollup/$BATCH     \
        $WORKING_DIR/UMI-count-rollup
done

#-------------------------------------------------------------------------------

cat UMI-count-rollup/* | gzip >prelim.UMI.count.rollup.gz

touch samples.to.exclude

echo "UMI_PER_CELL_CUTOFF = $UMI_PER_CELL_CUTOFF"

gunzip <prelim.UMI.count.rollup.gz \
| $DATAMASH_PATH -g 1 sum 3 \
| tr '|' '\t' \
| awk -v CUTOFF=$UMI_PER_CELL_CUTOFF '
    ARGIND == 1 {
        exclude[$1] = 1;
    } $3 >= CUTOFF && !($1 in exclude) {
        print $2 "\t" $1;
    }' samples.to.exclude - \
| sort -k1,1 -S 4G \
>final-output/cell.annotations

#rm final-output/UMI.count.matrix

gunzip < prelim.UMI.count.rollup.gz \
| tr '|' '\t' \
| awk '{
    if (ARGIND == 1) {
        gene_idx[$1] = FNR;
    } else if (ARGIND == 2) {
        cell_idx[$1] = FNR;
    } else if ($2 in cell_idx) {
        printf "%d\t%d\t%d\n",
            gene_idx[$3], cell_idx[$2], $4;
    }
}' final-output/gene.annotations final-output/cell.annotations - \
>final-output/UMI.count.matrix
#gawk needed

sudo chmod 777 /data/scratch/*


# # Rename files of fastq if needed before starting
# cd ./fastq
# for file in EB*; do
#     if [[ -f "$file" ]]; then
#         # Extract relevant parts using parameter expansion
#         prefix="${file%%_*}"
#         #echo $prefix
#         prefix="${prefix: -3}"
#         #echo $prefix
#
#         # Determine if the file contains "R1" or "R2"
#         if [[ $file == *"R1"* ]]; then
#             # Construct the new filename for R1
#             new_name="${prefix}_${SAMPLE_NAME}_R1_001.fastq.gz"
#             echo $new_name
#         elif [[ $file == *"R2"* ]]; then
#             # Construct the new filename for R2
#             new_name="${prefix}_${SAMPLE_NAME}_R2_001.fastq.gz"
#             echo $new_name
#         else
#             # If neither "R1" nor "R2" is found, skip the file
#             echo "not R1 or R2!"
#             continue
#         fi
#
#         # Rename the file
#         mv "$file" "$new_name"
#     fi
# done

#Rscript /home/endpipeline_sci.R
