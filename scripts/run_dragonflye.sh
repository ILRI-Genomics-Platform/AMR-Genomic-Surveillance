#!/bin/bash

#SBATCH -w compute06
#SBATCH --job-name=denovo_assembly
#SBATCH --output=stdout_%j.out
#SBATCH --error=stderr_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=batch


# load the required module(s)
module load dragonflye/1.2.1

# run dragonflye for denovo genome denovo_assembly


# define I/O

WORKDIR="/var/scratch/$USER/ACDC_AMR2025"
OUTDIR="${WORKDIR}/output/dragonflye"
fastq="${WORKDIR}/results/ont/klebsiella/nanoq/SRR28370682_filt.fastq.gz"

# check if output directory exists and create if not
if [ -d ${OUTDIR} ]; then
  echo -e "${OUTDIR}"
else
  mkdir -p ${OUTDIR}
fi


# run the assembly process
echo -e dragonflye \
    --reads ${fastq} \
    --gsize 5700000 \
    --prefix SRR28370682 \
    --outdir ${OUTDIR} \
    --assembler flye \
    --tmpdir ${WORKDIR} \
    --polypolish 1 \
    --minlen 500 \
    --mincov 2 \
    --force \
    --keepfiles \
    --depth 0 \
    --minreadlen 0 \
    --minquality 0 \
    --racon 1 \
    --medaka 0 \
    --namefmt "SRR28370682_%05d" \
    --cpus 4 \
    --ram 7 \
    --noreorient