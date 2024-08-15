#!/bin/bash

# Funzione per controllare se un file esiste
file_exists() {
    if [ -f "$1" ]; then
        return 0
    else
        return 1
    fi
}

# Scarica il genoma di riferimento e costruisci gli index se necessario
cd /genome/

# Autodetect il nome del file FASTA (con estensione .fa o .fasta)
genome_file=$(ls *.fa* | head -n 1)
filename=$(basename -- "$genome_file")
filename="${filename%.*}"

# Controlla se l'index per BWA esiste
bwa_index_check="/genome/bwa_${filename}/${filename}.fa.bwt"
if file_exists "${bwa_index_check}"; then
    echo "BWA index found. Skipping index generation."
else
    echo "BWA index not found. Generating index..."
    samtools faidx $genome_file

    # Build BWA MEM index
    mkdir bwa_$filename
    cd bwa_$filename
    cp ../$genome_file* ./
    bwa index $genome_file
    cd ..
fi

# Controlla se l'index per HISAT3N esiste
hisat_index_check="/genome/hisat3n_${filename}_CT/${filename}.fa.3n.CT.1.ht2"
if file_exists "${hisat_index_check}"; then
    echo "HISAT3N index found. Skipping index generation."
else
    echo "HISAT3N index not found. Generating index..."
    mkdir hisat3n_${filename}_CT
    cp ./$genome_file* ./hisat3n_${filename}_CT
    cd hisat3n_${filename}_CT
    /home/hisat-3n/hisat-3n-build --base-change C,T $genome_file $genome_file > hisat3n_${filename}_CT_index.log
    cd ..
fi

# Inizio della seconda parte dello script
threshold=$1

# Usa adattatori di default se non specificati
adapt1=$2
adapt2=$3
raw_dir="/scratch/raw.fastq"
# Trova tutti i file R1 con estensione fastq.gz o fq.gz
for in_fq_R1 in ${raw_dir}/*_R1.fastq.gz ${raw_dir}/*_R1.fq.gz ${raw_dir}/*_1.fastq.gz ${raw_dir}/*_1.fq.gz; do
    if [[ -f $in_fq_R1 ]]; then
        # Rileva l'estensione
        if [[ $in_fq_R1 == *.fastq.gz ]]; then
            ext="fastq.gz"
            base_name=$(basename "$in_fq_R1" .fastq.gz)
        elif [[ $in_fq_R1 == *.fq.gz ]]; then
            ext="fq.gz"
            base_name=$(basename "$in_fq_R1" .fq.gz)
        fi

        # Rimuovi il suffisso _R1, _R2, _1, _2 dal nome di base
        base_name=${base_name%_R1}
        sample=${base_name%_R2}

    # Definisci i percorsi per i file di input e output
    in_fq_R1=${raw_dir}/${sample}_R1.${ext}
    in_fq_R2=${raw_dir}/${sample}_R2.${ext}
    out_fq_R1=/scratch/fix.fastq/${sample}_R1_cutadapt.fastq.gz
    out_fq_R2=/scratch/fix.fastq/${sample}_R2_cutadapt.fastq.gz
    log=/scratch/fix.fastq/${sample}_cutadapt.log

    # Crea la directory per i file elaborati
    mkdir -p /scratch/fix.fastq
    cutadapt -j 0 --times 1 -e 0.1 -O 3 --quality-cutoff 25 -m 55 -a $adapt1 -A $adapt2 -o ${out_fq_R1} -p ${out_fq_R2} ${in_fq_R1} ${in_fq_R2} > ${log}

    # Esegui HISAT3N
    mkdir -p /scratch/bam.hisat3n
    in_fq_R1=/scratch/fix.fastq/${sample}_R1_cutadapt.fastq.gz
    in_fq_R2=/scratch/fix.fastq/${sample}_R2_cutadapt.fastq.gz
    out_bam=/scratch/bam.hisat3n/${sample}_hisat3n.bam
    ummapped_fq=/scratch/bam.hisat3n/${sample}_hisat3n_unmapped.fastq.gz
    log=/scratch/bam.hisat3n/${sample}_hisat3n.log
    ref_idx=/genome/hisat3n_${filename}_CT/$genome_file
    /home/hisat-3n/hisat-3n -x ${ref_idx} -1 ${in_fq_R1} -2 ${in_fq_R2} -p 20 --sensitive --base-change C,T --unique-only --repeat-limit 1000 --no-spliced-alignment -X 700 --un-conc-gz ${ummapped_fq} --summary-file ${log} --rg-id ${sample} --rg "PL:ILLUMINA" --rg "ID:"${sample} --rg "SM:"${sample} | samtools view -hb > ${out_bam}

    # Seleziona le letture con bassa qualità di mapping
    in_bam=/scratch/bam.hisat3n/${sample}_hisat3n.bam
    out_bam=/scratch/bam.hisat3n/${sample}_hisat3n.MAPQ20.LowerMAPQ20.bam
    samtools view -h -@ 4 ${in_bam} | awk '$1~"@" || $5 <= 20  {print $0}' | samtools view -@ 4 -hb > ${out_bam}

    # Ordina il BAM per nome di lettura
    in_bam=/scratch/bam.hisat3n/${sample}_hisat3n.MAPQ20.LowerMAPQ20.bam
    out_bam=/scratch/bam.hisat3n/${sample}_hisat3n.MAPQ20.LowerMAPQ20.SortName.bam
    temp_file=/scratch/bam.hisat3n/${sample}_hisat3n.MAPQ20.LowerMAPQ20.SortName.bam.temp
    samtools sort -O BAM -o ${out_bam} -T ${temp_file} -@ 15 -m 2G -n ${in_bam}

    # Estrai le letture con bassa qualità di mapping dal BAM
    in_bam=/scratch/bam.hisat3n/${sample}_hisat3n.MAPQ20.LowerMAPQ20.SortName.bam
    ref_genome_fa=/genome/hisat3n_${filename}_CT/$genome_file
    out_fq_R1=/scratch/bam.hisat3n/${sample}_hisat3n.LowerMAPQ20_R1.fastq.gz
    out_fq_R2=/scratch/bam.hisat3n/${sample}_hisat3n.LowerMAPQ20_R2.fastq.gz
    samtools fastq -@ 15 -0 /dev/null -s /dev/null -n -F 0x900 -1 ${out_fq_R1} -2 ${out_fq_R2} --reference ${ref_genome_fa} ${in_bam}

    # Unisci le letture non mappate e quelle con bassa qualità di mapping
    low_fq_R1=/scratch/bam.hisat3n/${sample}_hisat3n.LowerMAPQ20_R1.fastq.gz
    low_fq_R2=/scratch/bam.hisat3n/${sample}_hisat3n.LowerMAPQ20_R2.fastq.gz
    unmapped_fq_R1=/scratch/bam.hisat3n/${sample}_hisat3n_unmapped.fastq.1.gz
    unmapped_fq_R2=/scratch/bam.hisat3n/${sample}_hisat3n_unmapped.fastq.2.gz
    out_fq_R1=/scratch/bam.hisat3n/${sample}_R1_unmapped_and_LowerMAPQ20.fastq.gz
    out_fq_R2=/scratch/bam.hisat3n/${sample}_R2_unmapped_and_LowerMAPQ20.fastq.gz
    cat ${low_fq_R1} ${unmapped_fq_R1} > ${out_fq_R1}
    cat ${low_fq_R2} ${unmapped_fq_R2} > ${out_fq_R2}

    # Riallineamento con BWA MEM
    in_fq_R1=/scratch/bam.hisat3n/${sample}_R1_unmapped_and_LowerMAPQ20.fastq.gz
    in_fq_R2=/scratch/bam.hisat3n/${sample}_R2_unmapped_and_LowerMAPQ20.fastq.gz
    bwa_index=/genome/bwa_${filename}/$genome_file
    out_bam=/scratch/bam.hisat3n/${sample}_bwa_realign.bam
    bwa_log=/scratch/bam.hisat3n/${sample}_bwa_realign.log
    bwa mem ${bwa_index} ${in_fq_R1} ${in_fq_R2} -t 20 -M -R '@RG\tID:'${sample}'\tPL:ILLUMINA\tSM:'${sample} 2>${bwa_log} | samtools view -h -b -q 20 -f 3 -F 256 > ${out_bam}

    # Unisci i file BAM di HISAT-3n e BWA MEM
    in_bam_bwa=/scratch/bam.hisat3n/${sample}_bwa_realign.bam
    in_bam_hisat3n=/scratch/bam.hisat3n/${sample}_hisat3n.bam
    out_bam=/scratch/bam.hisat3n/${sample}_merge.MAPQ20.bam
    samtools cat -o ${out_bam} ${in_bam_hisat3n} ${in_bam_bwa}

    # Ordina il BAM per coordinate genomiche
    in_bam=/scratch/bam.hisat3n/${sample}_merge.MAPQ20.bam
    out_bam=/scratch/bam.hisat3n/${sample}_merge_sort.MAPQ20.bam
    temp_file=/scratch/bam.hisat3n/${sample}_merge_sort.MAPQ20.bam.temp
    samtools sort -O BAM -o ${out_bam} -T ${temp_file} -@ 15 -m 2G ${in_bam}

    # Rimuovi duplicati
    in_bam=/scratch/bam.hisat3n/${sample}_merge_sort.MAPQ20.bam
    out_log=/scratch/bam.hisat3n/${sample}_merge_sort_rmdup.log
    out_bam=/scratch/bam.hisat3n/${sample}_merge_sort_rmdup.MAPQ20.WithClip.bam
    out_matrix=/scratch/bam.hisat3n/${sample}_merge_sort_rmdup.matrix
    java -Xms50g -Xmx50g -XX:ParallelGCThreads=20 -jar /home/picard.jar MarkDuplicates I=${in_bam} O=${out_bam} M=${out_matrix} ASO=coordinate REMOVE_DUPLICATES=true 2> ${out_log}

    # Filtra clip, letture non concordanti, letture a bassa qualità MAPQ e allineamenti secondari
    in_bam=/scratch/bam.hisat3n/${sample}_merge_sort_rmdup.MAPQ20.WithClip.bam
    out_bam=/scratch/bam.hisat3n/${sample}_merge_sort_rmdup.MAPQ20.bam
    ref_genome_fa=/genome/hisat3n_${filename}_CT/$genome_file
    samtools view -@ 4 -h ${in_bam} -q 20 -f 3 -F 256 | /home/samclip --ref ${ref_genome_fa} --max 3 --progress 0 | awk 'function abs(v) {return v < 0 ? -v : v} $1~"@" || ($7 == "=" && abs($9) <= 2500 ) {print $0}' | samtools view -@ 4 -hb > ${out_bam}

    # Crea l'indice del file BAM
    in_bam=/scratch/bam.hisat3n/${sample}_merge_sort_rmdup.MAPQ20.bam
    out_bam_index=/scratch/bam.hisat3n/${sample}_merge_sort_rmdup.MAPQ20.bam.bai
    samtools index -@ 10 ${in_bam} ${out_bam_index}

    # Converte il file BAM in formato pmat
    mkdir -p /scratch/pmat_and_mpmat
    in_bam=/scratch/bam.hisat3n/${sample}_merge_sort_rmdup.MAPQ20.bam
    ref_genome_fa=/genome/hisat3n_${filename}_CT/$genome_file
    out_pmat=/scratch/pmat_and_mpmat/${sample}_CLEAN.pmat
    out_log=/scratch/pmat_and_mpmat/${sample}_CLEAN.log
    /root/anaconda3/envs/DetectSeq/bin/python /home/Detect-seq/src/detect_seq/bam2pmat.py -i ${in_bam} -r ${ref_genome_fa} -o ${out_pmat} -p 20 --out_format pmat --bed_like_format True --mut_type ALL --block_size 100000  --cover_num_cutoff 0 --mut_num_cutoff 0 --mut_ratio_cutoff 0 --keep_temp_file False --out_header False > ${out_log}

    # Definisci il locus di interesse
    chr="chr12"
    start=114437957
    end=114697959

    dir=$(dirname "$out_pmat")
    new_file="${sample}_IGH.pmat"
    new_path="${dir}/${new_file}"

    # Cerca nel file tutte le righe che iniziano con il locus di interesse
    grep "^${chr}\s" ${out_pmat} | awk -v chr="$chr" -v start="$start" -v end="$end" '$2>=start && $3<=end'  > ${new_path}

    directoryOut="/scratch/Res_Stat/${sample}"
    input_file=/scratch/pmat_and_mpmat/${sample}_IGH.pmat
    filenameL=$(basename "$input_file")
    filebase="${filenameL%.*}"
    path_without_filename=$(dirname "$input_file")
    search_strings=("CT" "GA")

    for search_string in "${search_strings[@]}"; do
        mkdir -p "${directoryOut}/filtered_${search_string}"
        awk -v dirOut="${directoryOut}/filtered_${search_string}" -v fname="$filenameL" -v threshold="$threshold" -v search="$search_string" '
            $9 ~ search {
                line = sprintf("%s_%s", $1, $2)
                for (i = 5; i <= NF; i++) {
                    if (i != 10 && i != 11 && i != 15 && i != 16) {
                        line = sprintf("%s %s", line, $i)
                    }
                }
                echo $13
                if ($13 >= threshold) {
                    file2 = dirOut "/filtered2_" fname
                    print line >> file2
                    close(file2)
                }
            }
        ' "$input_file" > "${directoryOut}/filtered_${search_string}/filtered_${filenameL}"

        awk -v search="$search_string" '
            {gsub("_", " "); if ($7=="CT") {print "variableStep chrom="$1" span=1\n"$2" "$9} else {print "variableStep chrom="$1" span=1\n"$2" -"$9}}
        ' "${directoryOut}/filtered_${search_string}/filtered_${filenameL}" >> "${directoryOut}/filtered_${search_string}/filtered_${filebase}.wig"

        awk -v search="$search_string" '
            {gsub("_", " "); if ($7=="CT") {print "variableStep chrom="$1" span=1\n"$2" "$9} else {print "variableStep chrom="$1" span=1\n"$2" -"$9}}
        ' "${directoryOut}/filtered_${search_string}/filtered2_${filenameL}" >> "${directoryOut}/filtered_${search_string}/filtered2_${filebase}.wig"

        output_file="${directoryOut}/filtered_${search_string}/filtered2_${filebase}.bed"

        while IFS= read -r line; do
            chr=$(echo "$line" | awk '{split($0, a, "[ _]"); print a[1]}')
            start=$(echo "$line" | awk '{split($0, a, "[ _]"); print a[2]}')
            count=$(echo "$line" | awk '{split($0, a, " "); print a[length(a)-1]}')
            for ((i = 0; i < count; i++)); do
                uniqueName=$i
                randomNumber=$(shuf -i 30-60 -n 1)
                strand=$(awk 'BEGIN{srand(); r = int(rand() * 2); if(r == 0) print "+"; else print "-"}')
                echo -e "$chr\t$start\t$start\t$uniqueName\t$randomNumber\t$strand" >> "$output_file"
            done
        done < "${directoryOut}/filtered_${search_string}/filtered2_${filenameL}"
    done

    for search_string in "${search_strings[@]}"; do
        echo "Processing $search_string"
        awk -f /home/locationAWK_COUNT.awk 114498650 114503581 12 SA "${directoryOut}/filtered_${search_string}/filtered_${filenameL}"
        awk -f /home/locationAWK_COUNT.awk 114498650 114503581 12 SA "${directoryOut}/filtered_${search_string}/filtered2_${filenameL}"
        awk -f /home/locationAWK_COUNT.awk 114661796 114665162 12 SU "${directoryOut}/filtered_${search_string}/filtered_${filenameL}"
        awk -f /home/locationAWK_COUNT.awk 114661796 114665162 12 SU "${directoryOut}/filtered_${search_string}/filtered2_${filenameL}"
        awk -f /home/locationAWK_COUNT.awk 114589303 114616242 12 SY3 "${directoryOut}/filtered_${search_string}/filtered_${filenameL}"
        awk -f /home/locationAWK_COUNT.awk 114589303 114616242 12 SY3 "${directoryOut}/filtered_${search_string}/filtered2_${filenameL}"
        awk -f /home/locationAWK_COUNT.awk 114568421 114581890 12 SY1 "${directoryOut}/filtered_${search_string}/filtered_${filenameL}"
        awk -f /home/locationAWK_COUNT.awk 114568421 114581890 12 SY1 "${directoryOut}/filtered_${search_string}/filtered2_${filenameL}"
        awk -f /home/locationAWK_COUNT.awk 114544419 114557888 12 SY2b "${directoryOut}/filtered_${search_string}/filtered_${filenameL}"
        awk -f /home/locationAWK_COUNT.awk 114544419 114557888 12 SY2b "${directoryOut}/filtered_${search_string}/filtered2_${filenameL}"
    done
fi
done
