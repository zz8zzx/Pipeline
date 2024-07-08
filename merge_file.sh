#!/bin/bash
#SBATCH --job-name=transcript_alignment
#SBATCH --output=/home/whuang/logs/%x_%j.log
#SBATCH --error=/home/whuang/logs/%x_%j.err
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

##############################################
# Navigate to project directory
cd "$project"

# Create the output directory for merged assemblies if it doesn't exist
mkdir -p merged_assemblies

# Loop through the merge lists
for merge_list in merge_lists/stringtie_merge_*.txt; do
    # Extract the group name from the merge list filename
    group=$(basename "$merge_list" | cut -d '_' -f 3 | cut -d '.' -f 1)

    # Create an array to hold the paths to the GTF files
    gtf_files=()

    # Read each sample name from the merge list and add its GTF path to the array
    while IFS= read -r sample_name; do
        gtf_files+=("stringtie_assemblies/${sample_name}.gtf")
    done < "$merge_list"

    # Run StringTie-Merge for the current group
    stringtie --merge -G ref/gencode.v44.basic.annotation.gtf -o "merged_assemblies/merged_${group}_transcripts.gtf" "${gtf_files[@]}"
done

# Loop through the merge lists
for merge_list in merge_lists/stringtie_merge_*.txt; do
    # Extract the group name from the merge list filename
    group=$(basename "$merge_list" | cut -d '_' -f 3 | cut -d '.' -f 1)

    # Sort the merged assembly GTF file
    igvtools sort "merged_assemblies/merged_${group}_transcripts.gtf" "merged_assemblies/merged_${group}_transcripts.sorted.gtf"

    # Index the sorted GTF file
    igvtools index "merged_assemblies/merged_${group}_transcripts.sorted.gtf"
done

################################################
## Check the exit status
if [ $? -ne 0 ]; then
  echo "merge_file failed at `date`."
  exit 1
else
  echo "merge_file completed successfully at `date`."
fi
