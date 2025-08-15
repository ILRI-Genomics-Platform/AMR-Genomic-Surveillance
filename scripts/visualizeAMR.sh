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
module load R/4.4
export R_LIBS=/export/apps/R/4.4/training-libs


# I/O
# WORKDIR="/var/scratch/$USER/ACDC_AMR2025"
# OUTDIR="${WORKDIR}/results/ont/klebsiella"

# Rscript ~/scripts/visualizeAMR.R \
#     --tree ${OUTDIR}/iqtree/core-snp.treefile \
#     --mlst ${OUTDIR}/mlst/*.tsv \
#     --pointfinder ${OUTDIR}/resfinder/*/PointFinder_results.txt \
#     --resfinder ${OUTDIR}/resfinder/*/ResFinder_results_tab.txt \
#     --prefix phylogeny-amr \
#     --outdir ${OUTDIR}/plots

Rscript ~/scripts/visualizeAMR.R \
    --tree ~/iqtree/core-snp.treefile \
    --mlst ~/mlst/*.tsv \
    --pointfinder ~/resfinder/*/PointFinder_results.txt \
    --resfinder ~/resfinder/*/ResFinder_results_tab.txt \
    --prefix phylogeny-amr \
    --outdir ~/plots
