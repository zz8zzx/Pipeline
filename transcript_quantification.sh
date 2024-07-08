#!/bin/bash
#SBATCH --job-name=transcript_quantification
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=30G
#SBATCH --time=04:00:00
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

# Reference url
reference_transcriptome_url=$2

# Create reference directory
ref_dir="$project/ref/"
mkdir -p "$ref_dir"

# Directory containing trimmed FASTQ files
trimmed_fastq_dir="$project/trimmed_fastq/"
quant_output_dir="$project/salmon_quant/"

# Create the output directory if it doesn't exist
mkdir -p "$quant_output_dir"

############################################################
# Navigate to project directory
cd "$project"

# Download reference transcriptome
cd "$ref_dir"
# curl -O "$reference_transcriptome_url"
reference_basename=$(basename "$reference_transcriptome_url")
gunzip "$reference_basename"
cd ..

# Index the transcriptome using Salmon
salmon index -p 6 -t "${ref_dir}$reference_basename" -i salmon_index

# Run Salmon quantification
for file in "$trimmed_fastq_dir"*_1_trimmed.fastq; do
  base="${file%_1_trimmed.fastq}"
  fastq2="${base}_2_trimmed.fastq"
  filename=$(basename "$base")
  sample_output_dir="${quant_output_dir}${filename}/"
  salmon quant -i salmon_index -l A -1 "$file" -2 "$fastq2" -p 6 --output "$sample_output_dir"
done

#############################################
## Check the exit status
if [ $? -ne 0 ]; then
  echo "transcript_quantification failed at `date`."
  exit 1
else
  echo "transcript_quantification completed successfully at `date`."
fi