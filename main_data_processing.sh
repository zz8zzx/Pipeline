#!/bin/bash
#SBATCH --job-name=main_data_processing
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --partition=work

############################################################
## Define input and output directories
# Set project directory
project="/scratch/peb004/whuang/$1"

# Set reference genome URL
reference_genome="$2"
reference_transcriptome="$3"
reference_gtf="$4"

# Pipeline directory
pipeline_dir="/home/whuang/pipeline/"

# Create project directory if it doesn't exist
mkdir -p "$project"

# Create log directory within the project directory
log_dir="$project/logs"
mkdir -p "$log_dir"

############################################################
# Function to log the status of each step
log_status() {
    if [ $? -eq 0 ]; then
        echo "$1 completed successfully at `date`" >> $log_dir/main.log
    else
        echo "$1 failed at `date`" >> $log_dir/main.log
        exit 1
    fi
}

#Call subroutines
# echo "Data Acquisition started at `date`" >> $log_dir/main.log
# bash ${pipeline_dir}data_acquisition.sh "$project" &> $log_dir/data_acquisition.log
# log_status "Data Acquisition"

echo "Data QC started at `date`" >> $log_dir/main.log
bash ${pipeline_dir}data_qc.sh "$project" &> $log_dir/data_qc.log
log_status "Data QC"

# echo "Transcript Quantification started at `date`" >> $log_dir/main.log
# bash ${pipeline_dir}transcript_quantification.sh "$project" "$reference_transcriptome" &> $log_dir/transcript_quantification.log
# log_status "Transcript Quantification"

# echo "Transcript Alignment started at `date`" >> $log_dir/main.log
# bash ${pipeline_dir}transcript_alignment.sh "$project" "$reference_genome" "$reference_gtf" &> $log_dir/transcript_alignment.log
# log_status "Transcript Alignment"

echo "Transcript Alignment (STAR) started at `date`" >> $log_dir/main.log
bash ${pipeline_dir}transcript_alignment_STAR.sh "$project" "$reference_genome" "$reference_gtf" &> $log_dir/transcript_alignment_STAR.log
log_status "Transcript Alignment (STAR)"
