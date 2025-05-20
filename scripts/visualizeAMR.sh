#!/bin/bash

#SBATCH -w compute06
#SBATCH --job-name=plot-phylo-amr
#SBATCH --output=stdout_%j.out
#SBATCH --error=stderr_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=batch


# load the required module(s)
module load R/4.2


# I/O
WORKDIR="/var/scratch/$USER/ACDC_AMR2025"
OUTDIR="${WORKDIR}/results/ont/klebsiella"

Rscript ~/scripts/visualizeAMR.R \
    --tree ${OUTDIR}/iqtree/core-snp.treefile \
    --mlst ${OUTDIR}/mlst/*.tsv \
    --pointfinder ${OUTDIR}/resfinder/*/PointFinder_results.txt \
    --resfinder ${OUTDIR}/resfinder/*/ResFinder_results_tab.txt \
    --prefix phylogeny-amr \
    --outdir ${OUTDIR}/plots