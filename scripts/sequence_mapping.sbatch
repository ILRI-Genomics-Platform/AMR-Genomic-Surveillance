#!/bin/bash
#SBATCH --job-name=bwa_align
#SBATCH --output=bwa_align_%j.log
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# Load modules
module load samtools/1.17 bwa/0.7.19

# Set variables
WORKDIR="/var/scratch/$USER/"
REFERENCE="./test_data/reference/genome.fa"
SAMPLE_ID="sample1"
READ1="./test_data/fastq/${SAMPLE_ID}_R1.fastq.gz"
READ2="./test_data/fastq/${SAMPLE_ID}_R2.fastq.gz"
OUTPUT_DIR="./test_data/output"

# Change into work directory and Create output directory
cd $WORKDIR
mkdir -p $OUTPUT_DIR

# Start timing
start_time=$(date +%s)

# Indexing the reference 
echo "Indexing reference genome in readiness for alignment"
bwa index $REFERENCE

# Run alignment
echo "Starting alignment for $SAMPLE_ID"
bwa mem -t $SLURM_CPUS_PER_TASK $REFERENCE $READ1 $READ2 | \
    samtools sort -@ $SLURM_CPUS_PER_TASK -o $OUTPUT_DIR/${SAMPLE_ID}.bam -

# Index BAM file
samtools index $OUTPUT_DIR/${SAMPLE_ID}.bam

# Calculate statistics
samtools flagstat $OUTPUT_DIR/${SAMPLE_ID}.bam > $OUTPUT_DIR/${SAMPLE_ID}.flagstat

# End timing
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Job completed in $elapsed_time seconds"
