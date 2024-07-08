#!/bin/bash
#SBATCH --job-name=data_acquisition
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --time=36:00:00
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
# Create directory for fastq files and download data using fasterq-dump
fastq_dir="$project/fastq/"
mkdir -p $fastq_dir

# SRR_Acc_List file name
srr_list=/home/whuang/input/$(basename "$project")_list.txt

##############################################################################
# Navigate to submitted direction and copy SRR_Acc_List.txt (RNA list file) to scratch project folder
cd $SLURM_SUBMIT_DIR
echo "Slurm submitted direction is "`pwd`
cp /home/whuang/input/SRR_Acc_List.txt $project

# Navigate to project directory
cd "$project"

# Run fasterq-dump with specified options
fasterq-dump --split-files --mem 5GB $(cat $srr_list) --progress --outdir $fastq_dir

#############################################
## Check the exit status
if [ $? -ne 0 ]; then
  echo "data_acquisition failed at `date`."
  exit 1
else
  echo "data_acquisition completed successfully at `date`."
fi