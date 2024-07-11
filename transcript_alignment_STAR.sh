#!/bin/bash
#SBATCH --job-name=transcript_alignment_STAR
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=12:00:00
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

# Reference URLs
reference_genome_url=$2
reference_gtf_url=$3

# Directory containing trimmed FASTQ files
trimmed_fastq_dir="$project/trimmed_fastq/"

# STAR output directories
star_index_dir="$project/star_index/"
mkdir -p "$star_index_dir"
star_output_dir="$project/star_align/"
mkdir -p "$star_output_dir"

# Create reference directory
ref_dir="$project/ref/"
mkdir -p "$ref_dir"

# Create the output directory for featureCounts results
featurecounts_output_dir="$project/featurecounts_results/"
mkdir -p "$featurecounts_output_dir"

##############################################
# Navigate to project directory
cd "$project"

# Navigate to reference directory and download reference genome
cd "$ref_dir"
reference_genome_url_basename=iwgsc_refseqv2.1_assembly.fa.gz
# Uncomment the following lines if you need to download and unzip the files
# wget "$reference_genome_url"
# gunzip "$reference_genome_url_basename"
# wget "$reference_gtf_url"
reference_gtf_url_basename=$(basename "$reference_gtf_url")
# gunzip "$reference_gtf_url_basename"
cd ..

# Build the STAR index
# STAR \
#   --runThreadN 32 \
#   --runMode genomeGenerate \
#   --genomeDir "$star_index_dir" \
#   --genomeFastaFiles "${ref_dir}${reference_genome_url_basename%.gz}" \
#   --sjdbGTFfile "${ref_dir}${reference_gtf_url_basename%.gz}" \
#   --sjdbOverhang 100 \
#   --limitGenomeGenerateRAM 80000000000

# Run STAR alignment
# for file in "$trimmed_fastq_dir"*_1_trimmed.fastq; do
#   base="${file%_1_trimmed.fastq}"
#   fastq2="${base}_2_trimmed.fastq"
#   sample_name=$(basename "$base")
#   mkdir -p "${star_output_dir}$sample_name"
#   STAR \
#     --runThreadN 32 \
#     --genomeDir "$star_index_dir" \
#     --readFilesIn "$file" "$fastq2" \
#     --outFileNamePrefix "${star_output_dir}$sample_name/" \
#     --outSAMtype SAM

  # Convert SAM to BAM and sort
#   samtools view -@ 32 -bS "${star_output_dir}$sample_name/Aligned.out.sam" | \
#     samtools sort -@ 32 -o "${star_output_dir}$sample_name/Aligned.sortedByCoord.out.bam"

  # Index the BAM file with `-c` option to handle large files
#   samtools index -c "${star_output_dir}$sample_name/Aligned.sortedByCoord.out.bam"
# done

# Run featureCounts for gene expression quantification
featureCounts -T 32 \
  -p \
  -a "${ref_dir}${reference_gtf_url_basename%.gz}" \
  -o "${featurecounts_output_dir}gene_expression_counts.txt" \
  "${star_output_dir}"*/Aligned.sortedByCoord.out.bam

################################################
## Check the exit status
if [ $? -ne 0 ]; then
  echo "transcript alignment (STAR) failed at `date`."
  exit 1
else
  echo "transcript alignment (STAR) completed successfully at `date`."
fi