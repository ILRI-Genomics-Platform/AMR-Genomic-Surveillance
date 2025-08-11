#!/usr/bin/bash
#SBATCH --job-name=r-pkgs
#SBATCH --output=r-pkgs-install_%j.log
#SBATCH -n 2
#SBATCH --mem=8G
#SBATCH --qos=normal
#SBATCH -p batch

module load R/4.4

Rscript install-r-packages.R  

 
