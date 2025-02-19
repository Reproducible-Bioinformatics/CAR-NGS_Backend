current_dir=$(pwd)
./generateLibInfo.sh "$current_dir"/Data/HTGTS.xlsx HTGTS_mouse
./htgts.sh "$current_dir"/Data/ filtered_JT8572_S1_R1_001.fastq filtered_JT8572_S1_R2_001.fastq libseqInfo.txt libseqInfo2.txt output_dir HTGTS_mouse mm9
