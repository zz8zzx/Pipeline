#!/bin/bash
#SBATCH --job-name=transcript_alignment
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=06:00:00
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
reference_genome_url=$2
reference_gtf_url=$3

# Directory containing trimmed FASTQ files
trimmed_fastq_dir="$project/trimmed_fastq/"

# Hisat2 output files
hisat2_index_dir="$project/hisat2_index/"
mkdir -p "$hisat2_index_dir"
hisat2_output_dir="$project/hisat2_align/"
mkdir -p "$hisat2_output_dir"

# Create reference directory
ref_dir="$project/ref/"
mkdir -p "$ref_dir"

# Create the output directory for StringTie assemblies if it doesn't exist
stringtie_output_dir="$project/stringtie_assemblies/"
mkdir -p "$stringtie_output_dir"

##############################################
# Navigate to project directory
cd "$project"

# Navigate to reference directory and download reference genome
cd "$ref_dir"
# wget "$reference_genome_url"
reference_genome_url_basename=$(basename "$reference_genome_url")
# gunzip "$reference_genome_url_basename"
# wget "$reference_gtf_url"
reference_gtf_url_basename=$(basename "$reference_gtf_url")
# gunzip "$reference_gtf_url_basename"
cd ..

# Build the Hisat2 index
hisat2-build -p 16 "${ref_dir}${reference_genome_url_basename%.gz}" "${hisat2_index_dir}hisat_index"

# Run Hisat2 alignment
for file in "$trimmed_fastq_dir"*_1_trimmed.fastq; do
  base="${file%_1_trimmed.fastq}"
  fastq2="${base}_2_trimmed.fastq"
  sample_name=$(basename "$base")
  mkdir -p "${hisat2_output_dir}$sample_name"
  hisat2 --rna-strandness RF --novel-splicesite-outfile "${hisat2_output_dir}$sample_name/splicesite.txt" \
    -S "${hisat2_output_dir}$sample_name/accepted_hits.sam" -p 10 -x "${hisat2_index_dir}hisat_index" -1 "$file" -2 "$fastq2"

  # Convert SAM to BAM
  samtools view -@ 5 -bS -h "${hisat2_output_dir}$sample_name/accepted_hits.sam" > "${hisat2_output_dir}$sample_name/accepted_hits.bam"

  # Sort and index the BAM file
  samtools sort -@ 5 "${hisat2_output_dir}$sample_name/accepted_hits.bam" -o "${hisat2_output_dir}$sample_name/accepted_hits.sorted.bam"
  samtools index "${hisat2_output_dir}$sample_name/accepted_hits.sorted.bam"
done

# Run StringTie assembly
for sample_dir in "$hisat2_output_dir"*/; do
  sample_name=$(basename "$sample_dir")
  stringtie -G "${ref_dir}${reference_gtf_url_basename%.gz}" -o "${stringtie_output_dir}$sample_name.gtf" -l "$sample_name" "${sample_dir}accepted_hits.sorted.bam"
done


# Loop through the merge lists
# for merge_list in merge_lists/stringtie_merge_*.txt; do
    # Extract the group name from the merge list filename
#     group=$(basename "$merge_list" | cut -d '_' -f 3 | cut -d '.' -f 1)

    # Create an array to hold the paths to the GTF files
#     gtf_files=()

    # Read each sample name from the merge list and add its GTF path to the array
#     while IFS= read -r sample_name; do
#         gtf_files+=("stringtie_assemblies/${sample_name}.gtf")
#     done < "$merge_list"

    # Run StringTie-Merge for the current group
#     stringtie --merge -G ref/gencode.v44.basic.annotation.gtf -o "merged_assemblies/merged_${group}_transcripts.gtf" "${gtf_files[@]}"
# done


################################################
## Check the exit status
if [ $? -ne 0 ]; then
  echo "transcript_alignment failed at `date`."
  exit 1
else
  echo "transcript_alignment completed successfully at `date`."
fi