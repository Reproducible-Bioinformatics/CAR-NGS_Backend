#!/bin/bash

# Parameters
bamsave=$1
genome_dir="/genome/"
fastq_dir="/scratch/"
fixed_fastq_dir="/scratch/fix.fastq"
output_dir="/scratch/output"
cellranger_exe="cellranger"  # Specify the exact path to the Cell Ranger executable
rscript_path="/scratch/combine_matrices.R"  # Specify the correct path to the R script
log_file="${fixed_fastq_dir}/rename_log.txt"

# Create the output and fixed FASTQ directories if they don't exist
mkdir -p $output_dir
mkdir -p $fixed_fastq_dir

# Initialize the log file
echo "Log delle trasformazioni FASTQ" > "$log_file"
echo "Data: $(date)" >> "$log_file"
echo "-----------------------------" >> "$log_file"

# Function to extract and normalize the file extension
normalize_extension() {
  if [[ "$1" == *.fq.gz ]]; then
    echo ".fastq.gz"
  elif [[ "$1" == *.fastq.gz ]]; then
    echo ".fastq.gz"
  fi
}

# Function to check and correct file names
fix_fastq_name() {
  fastq_file="$1"
  ext=$(normalize_extension "$fastq_file")
  base_name=$(basename "$fastq_file" | sed "s/.fq.gz$//" | sed "s/.fastq.gz$//")

  # Extract the sample name
  sample_name=$(echo "$base_name" | awk -F'_' '{print $1}')

  # Check for "_S" followed by a number
  if [[ "$base_name" =~ _S[0-9]+ ]]; then
    s_index=$(echo "$base_name" | grep -o '_S[0-9]\+')
  else
    s_index="_S1"
  fi

  # Check for "_L" followed by three digits
  if [[ "$base_name" =~ _L[0-9]{3} ]]; then
    lane=$(echo "$base_name" | grep -o '_L[0-9]\{3\}')
  else
    lane="_L001"
  fi

  # Check for "_R" followed by 1 or 2, or "_1"/"_2"
  if [[ "$base_name" =~ _R[12] ]]; then
    read_type=$(echo "$base_name" | grep -o '_R[12]')
  elif [[ "$base_name" =~ _1$ ]]; then
    read_type="_R1"
  elif [[ "$base_name" =~ _2$ ]]; then
    read_type="_R2"
  else
    read_type="_R1"
  fi

  # Check for "_001" at the end
  if [[ "$base_name" =~ _[0-9]{3}$ ]]; then
    chunk=$(echo "$base_name" | grep -o '_[0-9]\{3\}$')
  else
    chunk="_001"
  fi

  # Construct the new file name
  new_base_name="${sample_name}${s_index}${lane}${read_type}${chunk}.fastq.gz"
  new_fastq_file="${fixed_fastq_dir}/${new_base_name}"

  # Check if the symlink already exists
  if [ -e "$new_fastq_file" ]; then
    echo "Warning: $new_fastq_file already exists, skipping..."
  else
    # Create the symlink with the corrected name
    ln -s "$fastq_file" "$new_fastq_file"
    echo "Created symlink from $fastq_file to $new_fastq_file"

    # Log the transformation
    echo "File: $fastq_file -> $new_fastq_file" >> "$log_file"
  fi
}

# Normalize all FASTQ files in the given directory
fastq_files=($(find "$fastq_dir" -name "*.fastq.gz" -o -name "*.fq.gz"))

for fastq_file in "${fastq_files[@]}"; do
  fix_fastq_name "$fastq_file"
done

# Find the FASTA genome file
genome_fasta=$(find "$genome_dir" -maxdepth 1 -name "*.fa" -o -name "*.fasta" | head -n 1)

if [ -z "$genome_fasta" ]; then
  echo "Error: No FASTA file found in directory $genome_dir"
  exit 1
fi

# Find the GTF file
gtf_file=$(find "$genome_dir" -name "*.gtf" -o -name "*.gff" | head -n 1)

if [ -z "$gtf_file" ]; then
  echo "Error: No GTF/GFF file found in directory $genome_dir"
  exit 1
fi

# Create the index directory name based on the FASTA file name
genome_base=$(basename "$genome_fasta" .fa)
genome_base=$(basename "$genome_base" .fasta)
cellranger_index_dir="${genome_dir}/${genome_base}_cellranger_index/star/Genome"
cellranger_index_dirReal="${genome_dir}/${genome_base}_cellranger_index/"

if [ ! -f "$cellranger_index_dir" ]; then
  echo "Cell Ranger index not found, generating index..."
  mkdir -p $cellranger_index_dir
  $cellranger_exe mkref --genome=custom_ref \
                        --fasta=$genome_fasta \
                        --genes=$gtf_file \
                        --ref-version=1.0.0 \
                        --nthreads=8 \
                        --output-dir=$cellranger_index_dir
  echo "Cell Ranger index generated."
else
  echo "Cell Ranger index already exists, skipping generation."
fi

# Initialize a variable to collect output directories
output_dirs=""

# Group the FASTQ files by sample and determine paired-end or single-end
declare -A samples
for fastq_file in $(find "$fixed_fastq_dir" -name "*.fastq.gz"); do
  # Extract the base of the file name before "_S1"
  base_name=$(basename "$fastq_file" | sed -E 's/_S[0-9]+_L[0-9]{3}.*//')

  # Add the file to the sample list
  samples["$base_name"]+="$fastq_file "
done

# Process each sample separately
for sample_name in "${!samples[@]}"; do
  echo "Processing sample: $sample_name"

  # Create a specific output directory for this sample
  sample_output_dir="${output_dir}/${sample_name}_output"
  mkdir -p $sample_output_dir

  # Run Cell Ranger count
  $cellranger_exe count --transcriptome=$cellranger_index_dirReal \
                        --fastqs=$fixed_fastq_dir \
                        --sample=$sample_name \
                        --id=${sample_name}_output \
                        --output-dir=$sample_output_dir \
                        --create-bam=$bamsave

  # Modify barcodes to include the sample name
  for matrix_type in "filtered_feature_bc_matrix" "raw_feature_bc_matrix"; do
    barcodes_file="${sample_output_dir}/outs/${matrix_type}/barcodes.tsv.gz"
    zcat $barcodes_file | sed "s/$/-${sample_name}/" | gzip > "${sample_output_dir}/outs/${matrix_type}/barcodes_modified.tsv.gz"
    echo "Barcodes for sample $sample_name in ${matrix_type} modified."
  done

  # Add the sample output directory to the list
  output_dirs+="${sample_output_dir} "

done

# Prepare and print the R script command
rscript_command="Rscript $rscript_path $output_dirs"
echo "The following R script command will be executed:"
echo $rscript_command

# Run the R script to combine matrices
eval $rscript_command

echo "Single-cell RNA-seq analysis pipeline completed for all samples."
echo "Il log delle trasformazioni è stato salvato in $log_file."
