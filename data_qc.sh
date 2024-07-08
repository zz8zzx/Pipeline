#!/bin/bash
#SBATCH --job-name=data_qc
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=16:00:00
#SBATCH --partition=work

#######################################################
# Source the Conda initialization script
source /home/whuang/miniconda3/etc/profile.d/conda.sh

# Activate the environment
conda activate rna_seq_demo

# Verify the activation
mamba env list

############################################################
## Define input and output directories
# Project directory
project=$1

# Directory containing fastq files
fastq_dir="$project/fastq/"

# Output directory for FastQC results
output_dir="$project/fastqc/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Output directory for trimmed FastQ files
output_trimmed_dir="$project/trimmed_fastq/"

# Create the output directory if it doesn't exist
mkdir -p "$output_trimmed_dir"

################################################
# Navigate to project directory
cd "$project"

# Run FastQC
for file in "$fastq_dir"*_1.fastq; do
  base="${file%_1.fastq}"
  fastq2="${base}_2.fastq"
  fastqc -t 8 -o "$output_dir" "$file" "$fastq2"
done

# Run MultiQC
multiqc "$output_dir" -o "$output_dir"

# Loop through _1.fastq files and trim adapters using fastp
for file in "$fastq_dir"*_1.fastq; do
  base="${file%_1.fastq}"
  fastq2="${base}_2.fastq"
  filename="$(basename "$base")"
  output1="${output_trimmed_dir}${filename}_1_trimmed.fastq"
  output2="${output_trimmed_dir}${filename}_2_trimmed.fastq"
  fastp -i "$file" -I "$fastq2" -o "$output1" -O "$output2" --detect_adapter_for_pe --thread 4
done

################################################
## Check the exit status
if [ $? -ne 0 ]; then
  echo "data_qc failed at `date`."
  exit 1
else
  echo "data_qc completed successfully at `date`."
fi