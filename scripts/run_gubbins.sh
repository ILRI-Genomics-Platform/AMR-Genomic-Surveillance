#!/bin/bash

#SBATCH -w compute06
#SBATCH --job-name=recombination-gubbins
#SBATCH --output=stdout_%j.out
#SBATCH --error=stderr_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=batch


# load the required module(s)
module load gubbins/3.4


# I/O
WORKDIR="/var/scratch/$USER/ACDC_AMR2025"
OUTDIR="${WORKDIR}/results/ont/klebsiella/gubbins"


run_gubbins.py \
    --threads 2 \
    --prefix $OUTDIR/core-snp \
    --iterations 5 \
    --min-snps 3 \
    --min-window-size 100 \
    --max-window-size 10000 \
    --filter-percentage 25.0 \
    ${WORKDIR}/results/ont/klebsiella/snippy-core/core-snp-clean.full.aln