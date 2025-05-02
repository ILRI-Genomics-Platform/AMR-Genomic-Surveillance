# Bioinformatics Analysis for Antimicrobial Resistance Genomic Surveillance  

---  

###### **_Trainers_**: [John Juma](https://github.com/ajodeh-juma), [Kennedy Mwangi](https://github.com/wanjauk), [Daniel Ouso](https://github.com/ousodaniel) & [Gilbert Kibet](https://github.com/kibet-gilbert)

---

- **Table of Contents**
  - [Introduction](#introduction)
  - [Scope of the tutorial](#scope-of-the-tutorial)
  - [Background](#background)
  - [pre-requisites](#pre-requisites)
  - [Set up](#set-up)
  - [Analysis preparations](#analysis-preprations)
    - [Logging into the HPC](#logging-into-the-hpc)
    - [Project organisation](#project-organisation)
    - [Data retrieval](#data-retrieval)
  - [Bioinformatics analysis](#bioinformatics-analysis)
    - [ONT](#ont)
    - [Illumina](#illumina)
      - [Step 1: Loading modules](#step-1-loading-modules)
      - [Step 2: Data quality assessment](#step-2-data-quality-assessment)
      - [Step 3: Quality and adapater filtering](#step-3-quality-and-adapater-filtering)
      - [Step 4: Reads alignment to reference genome](#step-4-reads-alignment-to-reference-genome)
      - [Step 5: Generate genome consensus](#step-5-generate-genome-consensus)
      - [Step 6: De novo genome assembly](#step-6-de-novo-genome-assembly)
      - [Step 7: Quality assessment of the assembled contigs](#step-7-quality-assessment-of-the-assembled-contigs)
      - [Step 8: Annotatio of the assembled contigs](#step-8-annotatio-of-the-assembled-contigs)
      - [Step 9: Identification of AMR genes](#step-9-identification-of-amr-genes)
      - [Step 10: Identification of the virulence factors](#step-10-identification-of-the-virulence-factors)
  - [References](#references)

## Introduction  
*Escherichia coli* are ubiquitous bacteria found in the environment including the gut of humans and other animals and consists of numerous types and strains. Most types of *E. coli* are normal inhabitants of the gut of animals and do not cause diseases. However, some of the strains of *E. coli* have acquired genes that enable them to cause disease. These strains are commonly associated with food poisoning leading to diarrhoea and are referred to as diarrheagenic *E. coli* (DEC). Transmission occurs primarily through contaminated food but can also occur via person-to-person transmission, animal contact and water. For example, Shiga toxin-producing *E. coli* (STEC) serotype O157:H7 causes bloody diarrhoea and has previosuly been responsible for outbreaks worldwide.   

*Klebsiella pneumoniae* -  

## Scope of the tutorial  
In this workshop we will tackle, hands-on, the basic principles employed to generate consensus genome sequences of *E. coli* and *K. pneumoniae* and identify the serotypes involved in outbreaks, anti-microbial resistance genes (AMRs) and possible virulence factors they may have. 

## Background  


## Pre-requisites  
This module will come after the Introductory Linux module and therefore assumes familiarity with basic Linux command-line use. It also assumes you have an account and are operating in the ILRI computing cluster from a local Linux environment. 

## Set up  
You will use your personal computers computers to log into the ILRI HPC cluster, which operates on a Linux-operating system. Since we will be working from the remote servers, you will not need special setup for your personal laptops. However, you will need to install a program that enables you to log into the HPC.

## Analysis preprations

### Logging into the HPC  
To sign in to HPC using the following command, you will have been assigned a ***username*** that looks like `Bio4InfoXX` and a ***password***.
1. Replace `<user_name>` in the command with the provided username and execute (click enter). 
2. Enter the password and execute. ***Note:*** The password will be typed in the background but will not be visible to you as you type.
```
ssh <user_name>@hpc.ilri.cgiar.org
```
There are two nodes to choose from: `compute05`  and `compute06`. If your username (`Bio4InfoXX`) ends with an ***Odd Number*** (1,3,5,7,9) use `compute05` and if it ends with n ***even number*** (2,4,6,8,0) use `compute06`. Now let us secure a four of CPUs in one of the HPC nodes.  
>Compute05
```
interactive -w compute05 -c 2 -J amr-surveillance -p batch
```
>Compute06
```
interactive -w compute06 -c 2 -J amr-surveillance -p batch
```

### Project organisation  
We will start by setting up the project directory structure and then conduct the analysis stepwise. To setup a well-structured project directory we need to create some directories to store our data and scripts. We will be conducting our a anlysis from a directory in the `scratch` space of the HPC.  

1. *Create a directory using your username in the scratch:*
>**Note**

>Once inside the `hpc`, all instances of ```$USER``` will be equivalent to the hpc username that you were assigned. Your username, by default, is stored in a variable called `USER`. By using it, you will not have to type-in your username, rather, your shell will automatically pick your username which is the value stored in the `USER` variable. The `$` (dollar) character-prefix to a variable name is used to call the value of that variable.

```
mkdir -p /var/scratch/$USER
cd /var/scratch/$USER
```
1. *Create project directories:*
> **Note:** 

> We create a project directory `amr-surveillance` to store all that pertains to this tutorial/project. Within `amr-surveillance` we created `data` and subdirectories to store our input data and results from different analysis steps. We create `scripts` directory to store scripts/code that we genenrate or need in the analysis.

```
mkdir -p AMR-Genomic-Surveillance/{data,scripts}
cd AMR-Genomic-Surveillance/
mkdir -p ./data/{database,fastq,fastqc,fastp,spades,quast,bowtie,samtools}
```


### Data retrieval  

## Bioinformatics analysis  

### ONT 

### Step 1: Data quality assessment

```
# Remove Adapters
porechop \
    --input SRR25570560.fastq.gz \
    --format fastq \
    --threads 2 > adapter-ont.fq

# Quality filter
nanoq \
    --min-len 1000 \
    --min-qual 0 \
    --input  SRR25570560.fastq.gz 1> filt-ont.fq

# Compress
pigz -p 2 -c -n subsample-ont.fq > results/SRR25570560.fastq.gz


# Quality steps before and after QC
gzip -cd SRR25570560.fastq.gz | fastq-scan -g 5249449 > results/summary/SRR25570560-original.json
gzip -cd results/SRR25570560.fastq.gz | fastq-scan -g 5249449 > results/summary/SRR25570560-final.json



# Nanoplot
# Nanopore Plots
mkdir results/summary/SRR25570560-original results/summary/SRR25570560-final
NanoPlot \
    --threads 2 \
    --fastq SRR25570560.fastq \
    --outdir results/summary/SRR25570560-original/ \
    --prefix SRR25570560-original_


cp results/summary/SRR25570560-original/SRR25570560-original_NanoPlot-report.html
cp results/summary/SRR25570560-original_NanoPlot-report.html

NanoPlot \
    --threads 2 \
    --fastq results/SRR25570560.fastq.gz \
    --outdir results/summary/SRR25570560-final/ \
    --prefix SRR25570560-final_

# Tools
    bbduk: 39.19
    fastp: 0.24.0
    fastqc: 0.12.1
    fastq-scan: 1.0.1
    lighter: 1.1.3
    nanoplot: 1.44.1
    nanoq: 0.10.0
    pigz: 2.8
    porechop: 0.2.4
    rasusa: 2.1.0


```

### Step 2: Genome assembly

```
# assembly

dragonflye \
    --reads SRR25570560.fastq.gz \
    --gsize 5249449 \
    --outdir results \
    --assembler flye --tmpdir /var/scratch/jjuma/outdir-ecoli-ont
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
    --namefmt "SRR25570560_%05d" \
    --cpus 2 \
    --ram 7 \
    --noreorient


# Tools
    any2fasta: 0.4.2
    assembly-scan: 1.0.0
    bwa: 0.7.19-r1273
    dragonflye: 1.2.1
    flash: 1.2.11
    flye: 2.9.5-b1801
    medaka: 2.0.1
    megahit: 1.2.9
    miniasm: 0.3-r179
    minimap2: 2.28-r1209
    nanoq: 0.10.0
    pigz: 2.8
    pilon: 1.24
    racon: 1.5.0
    rasusa: 2.1.0
    raven: 1.8.3
    samclip: 0.4.0
    samtools: 1.21
    shovill: 1.1.0
    shovill-se: 1.1.0
    skesa: 2.5.1
    spades.py: 4.1.0
    velvetg: 1.2.10
    velveth: 1.2.10
    unicycler: 0.5.1


# assembly

gzip -c -d SRR25570560.fna.gz > SRR25570560.fna

export PROKKA_DBDIR=$(echo "$(which prokka | sed "s=/prokka==")/../db")
env
mkdir tmp_prokka/
TMPDIR=tmp_prokka/ bactopia-prokka \
    --evalue 1e-09 --coverage 80 --centre Bactopia \
    --cpus 2 \
    --prefix SRR25570560 \
    --locustag SRR25570560 \
    --proteins proteins.faa \
    SRR25570560.fna

    
rm -rf tmp_prokka/

```
### Step 3: Database preparation

```
# Make blastdb of contigs, genes, proteins
mkdir blastdb
cat SRR25570560/SRR25570560.fna | \
    makeblastdb \
    -dbtype "nucl" \
    -title "Assembled contigs for SRR25570560" \
    -out blastdb/SRR25570560.fna

cat SRR25570560/SRR25570560.ffn | \ 
    makeblastdb \
    -dbtype "nucl" \
    -title "Predicted genes sequences for SRR25570560" \
    -out blastdb/SRR25570560.ffn


cat SRR25570560/SRR25570560.faa | \
    makeblastdb \
    -dbtype "prot" \
    -title "Predicted protein sequences for SRR25570560" \
    -out blastdb/SRR25570560.faa

# Compress
tar -cvf - blastdb/ | gzip -c > SRR25570560/SRR25570560-blastdb.tar.gz


gzip SRR25570560/*.gff
gzip SRR25570560/*.gbk
gzip SRR25570560/*.fna
gzip SRR25570560/*.faa
gzip SRR25570560/*.ffn
gzip SRR25570560/*.sqn
gzip SRR25570560/*.fsa
gzip SRR25570560/*.tbl

# Cleanup intermediate files
rm -rf SRR25570560.fna blastdb/ results/*.pdb results/*.pjs results/*.pot results/*.ptf results/*.pto




# Tools
    makeblastdb: 2.16.0+
    prokka: 1.14.6


```

### Step 4: AMR detection

```
# amrfinderplus
gzip -c -d SRR25570560.fna.gz > SRR25570560.fna
gzip -c -d SRR25570560.faa.gz > SRR25570560.faa
gzip -c -d SRR25570560.gff.gz > SRR25570560.gff

mkdir database
tar -xzf amrfinderplus.tar.gz -C database
AMRFINDER_DB=$(find database/ -name "AMR.LIB" | sed 's=AMR.LIB==')

AMRFINDER_DB=$(find amrfinderplus.tar.gz/ -name "AMR.LIB" | sed 's=AMR.LIB==')

# Full AMRFinderPlus search combining results
amrfinder \
    --nucleotide SRR25570560.fna \
    --protein SRR25570560.faa \
    --gff SRR25570560.gff \
    --annotation_format prokka \
    --plus --ident_min -1 --coverage_min 0.5 --translation_table 11 \
    --database $AMRFINDER_DB \
    --threads 2 \
    --name SRR25570560 > SRR25570560.tsv

# Tools
    amrfinderplus: 4.0.19
    amrfinderplus-database: 4.0.19


# MLST 

mkdir database
tar -xzf mlst.tar.gz -C database
MLST_DB=$(find database/ -name "mlst.fa" | sed 's=blast/mlst.fa==')

MLST_DB=$(find mlst.tar.gz/ -name "mlst.fa" | sed 's=blast/mlst.fa==')

mlst \
    --threads 2 \
    --blastdb $MLST_DB/blast/mlst.fa \
    --datadir $MLST_DB/pubmlst \
    --scheme ecoli_achtman_4 --minid 95 --mincov 10 --minscore 50 \
    SRR25570560.fna.gz \
    > SRR25570560.tsv


mlst: 2.23.0
mlst-database: 2.23.0-20240325


```

### Illumina  

### Step 1: Loading modules  
### Step 2: Data quality assessment 
### Step 3: Quality and adapater filtering  
### Step 4: Reads alignment to reference genome    
### Step 5: Generate genome consensus  
### Step 6: De novo genome assembly  
### Step 7: Quality assessment of the assembled contigs  
### Step 8: Annotatio of the assembled contigs   
### Step 9: Identification of AMR genes  
### Step 10: Identification of the virulence factors  

## References  
