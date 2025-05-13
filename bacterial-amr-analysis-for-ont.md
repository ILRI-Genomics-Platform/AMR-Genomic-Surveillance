# Bioinformatics Analysis for Antimicrobial Resistance Genomic Surveillance: Long Reads; _Klebsiella pneumoniae_  

---  

###### **_Trainers_**: [John Juma](https://github.com/ajodeh-juma), [Kennedy Mwangi](https://github.com/wanjauk), [Gilbert Kibet](https://github.com/kibet-gilbert) & [Ouso Daniel](www.linkedin.com/in/ousodaniel)

---

## Table of Contents
  - [Overview](#overview)
  - [Learning Objectives](#learning-objectives)
  - [Prerequisites](#background)
  - [prerequisites](#prerequisites)
  - [Scope of the Tutorial](#scope-of-the-tutorial)
  - [Set-up](#set-up)
      - [Workshop Environment](#workshop-environment)
      - [Logging into the HPC](#logging-into-the-hpc)
      - [Compute Node](compute-node)
      - [Project organisation](#project-organisation)
  - [Bioinformatics Analysis](#bioinformatics-analysis)
      - [About the Sample](#about-the-sample)
      - [Step 1: Data Quality Assessment](#step-1-data-quality-assessment)
      - [Step 2: Genome Assembly](#step-2-genome-assembly)
      - [Step 3: Genome Annotation](#step-3-genome-annotation)
      - [Step 4: AMR Detection](#step-4-amr-detection)
          - [Output Formart](#output-format)
          - [AMR Dectection with ResFinder](#amr-detection-with-resfinder)
          - [AMR Dectection with CARD/RGI](#amr-detection-with-card/rgi)
      - [step 5: Pathogen Relatedness](#step-5-pathogen-relatedness)
          - [MLST](#mlst)
              - [MLST Output Format](#mlst-output-format)
              - [MLST Results Interpretation](#mlst-results-interpretation)
              - [Visualising MLST Results](#visualising-mlst-results)
              - [MLST Interpretation Limitations](#mslt-interpretation-limitations)
  - [Step 6: Variant Calling and Consensus Assemblies](#step-6-variant-calling-and-consensus-assemblies)
      - [Fast Bacterial Variant Calling with Contigs](#fast-bacterial-variant-calling-with-contigs)
          - [Snippy Outputs](#snippy-outputs)
          - [Visualising Snippy Variants](#visualising-snippy-variants)
          - [Build Core and Whole Genome Alignments from Snippy Output](#build-core-and-whole-genome-alignments-from-snippy-output)
          - [Run Snippy-core](#run-snippy-core)
          - [Snippy Core Outputs](#snippy-core-outputs)
          - [Cleanup the Snippy SNP Alignment Intermediates](#cleanup-the-snippy-snp-alignment-intermediates)
          - [Compute Pairwise SNP Distances](#compute-pairwise-snp-distances)
- [Step 7: Dealing with Recombination](#step-7-dealing-with-recombination)
    - [Phylogenetic Analysis of Gubbins Output](#phylogenetic-analysis-of-gubbins-output)
    


## Overview
Welcome to this workshop on antimicrobial resistance (AMR) genomics analysis!

This workshop is designed for professionals who are already familiar with AMR and genomics concepts but are looking to develop practical bioinformatics skills for analysing genomic data related to AMR.
We'll explore the fundamental concepts, tools, and workflows used in AMR genomics analysis, with a focus on practical applications using _Escherichia coli_ and _Klebsiella pneumoniae_ sequencing data.

## Learning Objectives
By the end of this workshop, you will be able to:

- Understand the role of genomics in AMR surveillance and research
- Process and analyze logn-read sequencing data to identify AMR genes and mutations
- Use common bioinformatics tools to interpret AMR-related genomic data
- Implement basic genomic analysis workflows for AMR detection

## Prerequisites
This workshop assumes:

- Familiarity with basic microbiology and antimicrobial resistance mechanisms
- Understanding of fundamental genomics concepts
- Basic command-line skills

## Introduction  
<!-- *Escherichia coli* are ubiquitous bacteria found in the environment including the gut of humans and other animals and consists of numerous types and strains. Most types of *E. coli* are normal inhabitants of the gut of animals and do not cause diseases. However, some of the strains of *E. coli* have acquired genes or gen mutations that enable them to cause disease. These strains are commonly associated with food poisoning leading to diarrhoea and are referred to as diarrheagenic *E. coli* (DEC). Transmission occurs primarily through contaminated food but can also occur via person-to-person transmission, animal contact and water. For example, Shiga toxin-producing *E. coli* (STEC) serotype O157:H7 causes bloody diarrhoea and has previosuly been responsible for outbreaks worldwide.    -->

_Klebsiella pneumoniae_ is a Gram-negative, opportunistic pathogen notorious for causing hospital-acquired infections, including pneumonia, bloodstream infections, and urinary tract infections. A major concern is its propensity for antimicrobial resistance (AMR), particularly through extended-spectrum β-lactamases (ESBLs) and carbapenemases (e.g., KPC, NDM), which render last-resort antibiotics like carbapenems ineffective. Multidrug-resistant (MDR) and extensively drug-resistant (XDR) strains are increasingly reported worldwide, fueled by plasmid-mediated resistance genes and clonal spread in healthcare settings.

The rise of AMR in _K. pneumoniae_ underscores the urgent need for enhanced surveillance, infection control, and novel therapeutics. 

## Scope of the Tutorial  
In this workshop we will tackle, hands-on, the basic principles employed to generate consensus genome sequences of *E. coli* and *K. pneumoniae* and identify the serotypes involved in outbreaks, anti-microbial resistance genes (AMRs) and possible virulence factors they may have. 

## Set-up  
### Workshop Environmnet
The workshop is mainly based on working in a HPC envirnment. However, local set-ups will be provisioned, although some steps of the workflows may be impossible to execute on laptops. For access to the HPC, you will use your personal computers to log into the ILRI HPC cluster, which operates on a Linux operating system. If using a non-Unix operating system you will require a program that enables interfacing with the Linux HPC environment.

<!-- ## Analysis preprations -->

### Logging into the HPC  
To log in to the HPC, you use the provided details (username and password) with the below command. The username follows the parttern `Bio4InfoXX`, where `XX` is a number.
```
ssh <user_name>@hpc.ilri.cgiar.org
```
- Replace `<user_name>` in the command with the provided username and execute (press enter). 
- Next, you will be prompted for the password.  
 >***Note:*** The password will not be visible to you as you type--just a little faith.

### Compute Node
There are two nodes to choose from: `compute05`  and `compute06`. If your username (`Bio4InfoXX`) ends with an ***Odd Number*** (1,3,5,7,9...) use `compute05` and if it ends with an ***even number*** (2,4,6,8...) use `compute06`. For the purposes of this workshop we will each avail four vCPUs (two CPU cores).  
>Compute05
```
interactive -w compute05 -c 2 -J amr-surveillance -p batch
```
>Compute06
```
interactive -w compute06 -c 2 -J amr-surveillance -p batch
```
- `-w compute05`: specifically request to run the session on `compute05`.
- `-c 2`: allocate 2 CPU cores.
- `-J amr-surveillance`: names the job _amr-surveillance_.
- `-p batch`: specifies the `SLURM` partition (queue) to use, here `batch`.

<details>
    <summary>
        Click to toggle <b style='color:blue'>Modules to load</b>
    </summary>
</details>

### Project Organisation  
For any Bioinformatics project, it's good practice to have a structured file system to ensure proper logical separation of components. By so doing, we enhance clarity, collaboration and reproducibility. Before commencing our analysis, we will start by setting up the project directory structure. Widely adopted project structure include the following directories:
- `data`: stores copy of the original raw data
- `scripts`: stores utility scrippts and code
- `results`/`outputs`: stores analysis results/output
- `logs`: stores reference logs
- `tmp`/`scratch`: stores tempoary/intermediate process outputs

1. *Create a course directory called `ACDC_AMR2025`:*

```
# Creating the entity and project directory
mkdir -p /var/scratch/$USER/ACDC_AMR2025

# Change directory into the created directory
cd /var/scratch/$USER/ACDC_AMR2025
```
>**Note** Once inside the `hpc`, all instances of ```$USER``` will be equivalent to the `hpc` username that you were assigned. Your username, by default, is stored in a variable called `USER`. By using it, you will not have to type your username, rather, your shell will automatically retrieve your username which is the value stored in the `USER` variable. The `$` (dollar) prefix is used to retrieve the value of that variable.

2. *Create project sub-directories:*

```
# Create the project directory and sub-directories
mkdir -p AMR-Genomic-Surveillance/{data,scripts}

# Create tool sub-directories within the data directory
mkdir -p \
results/ont/klebsiella/{porechop,nanoq,fastq-scan,nanoplot,dragonflye,prokka,amrfinder,mlst,tmp/{dragonflye,prokka,amrfinder,snippy},snippy,snippy-core,gubbins}

ln -sf /var/scratch/global/jjuma/ACDC_AMR2025/[dpsr]* .
```
> **Note:**  We create a project root directory `ACDC_AMR2025` to store all that pertains to this tutorial/project. Within `ACDC_AMR2025` we created sub-directories aligned with the workflow tools.

## Bioinformatics Analysis

### About the Sample
**Project accession:** [PRJNA1087001](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1087001/)

**Project paper:** [Publication](https://elifesciences.org/reviewed-preprints/98300)

**Project summary:** Benchmark various nanopore-based variant callers on 14 different species. Samples are sequenced on the latest (September 2023) R10.4.1 Nanopore flowcells and Illumina. Ground truth assemblies are generated for each sample.

### Step 1: Data Quality Assessment
We start by exploring the quality of the raw sequence reads. We need to know whether the data is worth commencing analysis on, or otherwise.
```
# Initial quality assessment
NanoPlot
    --threads 2 # number of processing
    --fastq ./data/klebs/ont/SRR28370682.fastq.gz # input reads
    --outdir ./results/ont/klebsiella/nanoplot/ # output dir
    --prefix SRR28370682-original_ # output naming prefix
```
The HPC does not provide a graphical interface, therefore, we must transfer the quality assessment results to our local computer to view the outputs. The `rsync` or `scp` commands can help us with the transfer. From you local machine, and n the terminal, transfer the QC file as below.
```
# Transfer QC reports fro viewing
rsync -avzP \
    <hpc_user_name>@hpc.ilri.cgiar.org:<path_to_files_in_hpc> <local_path>
```
Sequencing adapters are short sequences added to genomic fragments during library preparation to enable sequencing, which must be removed before downstream analysis.


```
# Remove Adapters
porechop \
    --input SRR28370682.fastq.gz \ # reads
    --format fastq \ # read format
    --threads 2 > adapter-ont.fq # number of threads
```
We can filter ours reads based on the quality assessment observation.
```
# Quality filter
nanoq \
    --min-len 1000 \ # minimum read length
    --min-qual 0 \ # minimum read quality
    --input ./results/ont/klebsiella/porechop/SRR28370682_adapter_ont.fastq.gz \
    --output-type g \ # gzip compressed output
    --output ./results/ont/klebsiella/nanoq/SRR28370682_filt.fastq.gz
```


For a comprehensive quality visualization and reporting, we use `NanoPlot` again.
```
# Nanopore Plots
NanoPlot \
    --threads 2 \ # number of processing threads
    --fastq ./results/ont/klebsiella/nanoq/SRR28370682_filt.fastq.gz # input fastq file
    --outdir ./results/ont/klebsiella/nanoplot/ \ # output directory
    --prefix SRR28370682-filt_ # output name prefix
```
Let's copy the second round of quality assessment results to our local computers, as before.


### Step 2: Genome Assembly

How do we get genomes from reads?
Dragonflye is a pipeline for  quick and easy assembling of bacterial genomes from Oxford Nanopore reads.

<details>
  <summary>Click to toggle <span style="color:blue"><b>Dragonflye assembly steps</b></span></summary>

1. Estimate genome size and read length from reads (unless `--gsize` provided) (`kmc`)  
2. Filter reads by length (default `--minreadlength 1000`) (`nanoq`)  
3. Reduce FASTQ files to a sensible depth (default `--depth 150`) (`rasusa`)  
4. Remove adapters (requires `--trim` be given) (`porechop`)  
5. Assemble with `Flye`, `Miniasm`, or `Raven`  
6. Polish assembly with `Racon` and/or `Medaka`  
7. Polish assembly with short reads via `Polypolish` and/or `Pilon`  
8. Remove contigs that are too short, too low coverage, or pure homopolymers  
9. Produce final FASTA with nicer names and parsable annotations  
10. Reorient contigs from final FASTA using `dnaapler`  
11. Output parsable assembly statistics (`assembly-scan`)  

</details>

```
# Assembly

dragonflye \
    --reads ./results/ont/klebsiella/nanoq/SRR28370682_filt.fastq.gz \ # input reads
    --gsize 5248520 \ # reference genome size in bp
    --prefix SRR28370682 \
    --outdir ./results/ont/klebsiella/dragonflye \ # output directory
    --assembler flye # specify the assembly algorithm; other options miniasm, raven
    --tmpdir ./results/ont/klebsiella/tmp/dragonflye/
    --polypolish 1 \ # enable assembly polishing
    --racon 1 \ # enable racon consensus refinement
    --medaka 0 \ # disable medaka polishing
    --minlen 500 \ # min contig length in final assembly
    --mincov 2 \ # min contig coverage for retention
    --minreadlen 0 \ # disable filtering by length
    --minquality 0 \ # disable filtering by quality
    --force \ # overwrite existing out dir, if present
    --keepfiles \ # keep intermeadiates; for debugging
    --depth 0 \ # disable read sub-sampling; use all reads
    --noreorient # disable contig re-orientation based on reference    
    --namefmt "SRR25570560_%05d" \ # output formarting; how many digits
    --cpus 2 \ # number of processing cpus
    --ram 7 \ # amount of memory to reserve
```
    
We need to verify the quality of the resulting assembly before annotation. We will use the utility script below:
```
stats.sh in=./results/ont/klebsiella/dragonflye/SRR28370682.fa
```
   > Note Unfortunately, the N50 and L50 values generated by `stats.sh` are switched. N50 should be a length and L50 should be a count. The results table below shows the corrected values based on `stats.sh` outputs.

### Step 3: Genome Annotation

After assembling our genome, we need to characterise it through annotation to enable inferencing. Prokka (Prokaryotic annotation) is a pipeline tool designed for the rapid annotation of prokaryotic genomes, including bacteria and archaea. After Prokka annotation, tools like ABRicate or RGI can be run on the annotated genome to identify resistance genes in their proper genomic context.

Prokka respects the default Linux temporary directory--`TMPDIR`. Therefore we assign one if none exists.
```
# Only assigns /tmp to TMPDIR if it is unset
${TMPDIR:=/tmp}
```
However, because we are all running the same pipeline and we don't want to overwrite each other's files, we will use our project-specific temporary directories.
```
export TMPDIR=./results/ont/klebsiella/tmp/prokka/
```
The tool also comes with its own database which is stored in the variable `PROKKA_DBDIR`. To point it to a custom DB, the variable can be tailored and exported.
```
# To use a custom DB
export PROKKA_DBDIR=<path/to/custom/db>
```
Running prokka.
```
# Run prokka
prokka \
    --evalue 1e-09 # significance of alignment
    --coverage 80 # min alignment coverage
    --centre ILRI \ # metadata, centre name
    --cpus 2 \ # processing cpus
    --prefix SRR28370682 \ # output prefix
    --locustag SRR28370682 \ # gene tags
    --proteins ./databases/prokka/proteins.faa \ # custom priority to use for annotation befor PROKKA_DBDIR
    --outdir ./results/ont/klebsiella/prokka # output directory
    ./results/ont/klebsiella/dragonflye/SRR28370682.fa # input genome to be annotated

```

Prokka generates multiple output files in standard bioinformatics formats:

|File Extension | Format | Description|
|--- | --- |---|
`.gff` | GFF3 | General feature format (coordinates of all features)
`.gbk` | GenBank | Standard sequence with annotations format
`.fna` | FASTA | Nucleotide FASTA file of the input contig sequences
`.faa` | FASTA | Protein FASTA file of the translated CDS sequences
`.ffn` | FASTA | Nucleotide FASTA file of all feature sequences
`.sqn` | Sequin | ASN.1 format for GenBank submission
`.fsa` | FASTA | Nucleotide FASTA file of the input contig sequences with annotations
`.tbl` | Feature Table | Feature table for GenBank submission
`.txt` | Text | Statistics of the annotation run

We can clean intermediate files to conserve space:
```
rm -r ./results/ont/klebsiella/prokka/*.pdb \
./results/ont/klebsiella/prokka/*.pjs \
./results/ont/klebsiella/prokka/*.pot \
./results/ont/klebsiella/prokka/*.ptf \
./results/ont/klebsiella/prokka/*.pto

```

### Step 4: AMR Detection
Now that we have an annotated genome, we can query it for antimicrobial resistance. There a variety of tools for the task. We will use a leading and frequently updated one by NCBI--[***AMRFinderPlus***](https://github.com/ncbi/amr/wiki/Home). It identifies acquired antimicrobial resistance genes in bacterial protein and/or assembled nucleotide sequences as well as known resistance-associated point mutations for several taxa.

The most critical part of any AMR tool is the underlying database. So we take time to understand the databse on which AMRFinderPlus is based.

"This database is derived from the [Pathogen Detection Reference Gene Catalog](https://www.ncbi.nlm.nih.gov/pathogens/isolates#/refgene/), [Pathogen Detection Reference Gene Hierarchy](https://www.ncbi.nlm.nih.gov/pathogens/genehierarchy/), and [Reference HMM Catalog](https://www.ncbi.nlm.nih.gov/pathogens/hmm/) and is used by the [Pathogen Detection](https://ncbi.nlm.nih.gov/pathogens/) isolate analysis system to provide results to the [Isolates browser](https://www.ncbi.nlm.nih.gov/pathogens/isolates) and [MicroBIGG-E](https://www.ncbi.nlm.nih.gov/pathogens/microbigge) as well as the command-line version of AMRFinderPlus. The 'core' subset version focuses on acquired or intrinsic AMR gene products including point mutations in a limited set of taxa. As of version 4.0 AMRFinderPlus also includes [StxTyper](https://github.com/ncbi/stxtyper) which has a separate DNA-sequence database and algorithm for typing Stx operons.

The 'plus' subset include a less-selective set of genes of interest including genes involved in virulence, biocide, heat, metal, and acid resistance."

>**Note:** that AMRFinderPlus reports gene and point mutation presence/absence; it does not infer phenotypic resistance. Many of the resistance genes detected by AMRFinderPlus may not be relevant for clinical management or antimicrobial surveillance, AMRFinderPlus.

Let's get started by defining a temporary directory for the tool and define DB path:
```
export TMPDIR=./results/ont/klebsiella/tmp/amrfinder/

AMRFINDER_DB=./databases/amrfinderplus/2023-11-15.1/
```
You can update the database before the analysis to ensure the latest information.
```
# Updating the database
amrfinder --update
```
```
# Run AMR finder
amrfinder \
    --nucleotide ./results/ont/klebsiella/prokka/SRR28370682.fna \ # prokka out genomic
    --protein ./results/ont/klebsiella/prokka/SRR28370682.faa \ # prokka out proteins
    --gff ./results/ont/klebsiella/prokka/SRR28370682.gff \ # gene feature file
    --annotation_format prokka \ # specify input format producer
    --organism Klebsiella_pneumoniae \ # exploit prior knowledge to narrow the search
    --plus \ # detect/report beyond AMR genes
    --ident_min -1 \ # matc identity; any
    --coverage_min 0.5 \ # extent of gene to be covered
    --translation_table 11 \ # use bacterial codon table for DNA-Prot translation
    --database $AMRFINDER_DB \
    --threads 2 \
    --name SRR28370682 > ./results/ont/klebsiella/amrfinder/SRR28370682.tsv
```
#### Output Format
AMRFinderPlus produces a tab-separated output file with detailed information:
Column | Description
---|---
Protein identifier | Protein/contig identifier in the input sequence
Gene symbol | Gene symbol for the AMR gene
Sequence name | Matching sequence name in the AMR database
Scope | Core, Plus (e.g., virulence factors)
Element type | AMR gene class (e.g., AMR, STRESS, VIRULENCE)
Element subtype | Specific resistance mechanism
Class | Antibiotic class
Subclass | Specific antibiotic subclass
Method | Detection method used (PARTIALP, EXACTP, BLASTX, HMM)
Target length | Reference sequence length
Reference sequence length | Length of matching reference
% Coverage of reference | Percentage of reference covered by alignment
% Identity to reference | Percent identity to reference sequence
Alignment length | Length of the alignment
Accession of closest sequence | Accession number of the matching sequence
Name of closest sequence | Name of the matching sequence
HMM id | Identifier of HMM used for detection (if applicable)
HMM description | Description of HMM (if applicable)

#### AMR Detection with ResFinder
[**ResFinder**](https://genepi.dk/resfinder) identifies acquired genes and/or finds chromosomal mutations mediating antimicrobial resistance in total or partial DNA sequence of bacteria.
<details>
<summary>
    Click to toggle <b style='color:blue'>ResFinder CLI</b>
</summary>
We get started by availing the relevant databases (pre-downloaded).
```
# Creating the database directory
mkdir -p ./databases/resfinder/database
```

```
# Copy pre-available databases to local path
cp -rf
/var/scratch/global/gkibet/ilri-africa-cdc-training/wastewater-surveillance-analysis/database/{disinfinder_db,pointfinder_db,resfinder_db}
./databases/resfinder/database

# Export DB path variable
export DB_PATH_RES=./databases/resfinder/database/
```

    ```
    python -m resfinder \
    -ifa ./results/ont/klebsiella/dragonflye/SRR28370682.fa \ # input FASTA
    -o ./results/ont/klebsiella/resfinder/ \ # output
    -s klebsiella \ # spicies
    --min_cov 0.6 \ # min aligmnt gene covarage
    --threshold 0.9 \ # min id for amr match
    --min_cov_point 0.6 \ # min coverage for point mutations
    --threshold_point 0.9 \ # min id for point mutation
    --ignore_stop_codons \ # ignore premature stops; segmented assemblie
    --ignore_indels \ # inore insertion/deletions
    --acquired \ # acquired amr genes
    --point # find point mutations known for amr
    ```
</details>


<details>
    <summary>
        Click to toggle <b style='color:blue'>Batch ResFinder detection</b>
    </summary>

```
for fn in ./data/klebs/pathogenwatch/assemblies-to-test/*.fasta; do
    sample=$(basename $fn)
    sample="${sample%.*}"
    echo -e "-------------------------------\n"
    echo -e "running ResFinder on: $sample - $fn"
    
    python -m resfinder \
    -ifa $fn \ # input FASTA
    -o ./results/ont/klebsiella/resfinder/${sample} \ # output
    -s klebsiella \ # spicies
    --min_cov 0.6 \ # min aligmnt gene covarage
    --threshold 0.9 \ # min id for amr match
    --min_cov_point 0.6 \ # min coverage for point mutations
    --threshold_point 0.9 \ # min id for point mutation
    --ignore_stop_codons \ # ignore premature stops; segmented assemblie
    --ignore_indels \ # inore insertion/deletions
    --acquired \ # acquired amr genes
    --point # find point mutations known for amr
done
```
</details>

#### AMR Detection with CARD/RGI
[**CARD/RGI**](https://card.mcmaster.ca/analyze/rgi) can be used to predict resistomes from protein or nucleotide data based on homology and SNP models.

### Step 5: Pathogen Relatedness
Understanding the differences and relationships within circulating pathogens is an important aspect of genomics epidemiology. Whole-genome comparison has better resolution for pathogen characterisation than fragments per genome equivalent (FPGE) or gene-based or multi-locus sequence typing (MLST). Distances between genomes can be compared based a reference or _de novo_; _k-mer_-composition-based and _core-genome_ assembly-based.

#### MLST 

Multi-Locus Sequence Typing (MLST) is an essential tool in analyzing multiple antimicrobial resistant (MAR) bacterial isolates. It provides a standardized approach to characterizing bacterial strains based on the sequences of multiple housekeeping genes. It focuses on bacterial population structure and evolutionary relationships.

MLST is a molecular typing method that characterizes bacterial isolates by sequencing internal fragments (typically 450-500 bp) of multiple housekeeping genes (usually 7-8 genes). Each unique sequence for a gene is assigned an allele number, and the combination of allele numbers defines the sequence type (ST) of an isolate.

Create a database directory and extract a pre-downloaded MLST databse into it:
```
# Creating the database directory
mkdir -p ./databases/mlst/database
tar -xzf ./databases/mlst/mlst.tar.gz -C ./databases/mlst/database

# MLST base dir
MLST_DB_DIR=./databases/mlst/database/
```

```
mlst \
    --threads 2 \ # num of processing threads
    --blastdb $MLST_DB_DIR/blast/mlst.fa \ # specific BLAST-indexed database file containing all MLST allele sequences
    --datadir $MLST_DB_DIR/pubmlst \ # directory containing schem definitions
    --scheme klebsiella \ # specify scheme: gapA, infB, mdh, pgi, phoE, rpoB, and tonB for kleb
    --minid 95 \ # min match identity
    --mincov 10 \ # min match coverage
    --minscore 50 \ # min align score
    ./results/ont/klebsiella/prokka/SRR28370682.fna \ # input
    > ./results/ont/klebsiella/mlst/SRR28370682.tsv # output
```

<details>
    <summary>
        Click to toggle <b style='color:blue'>Batch MLST typing</b>
    </summary>

```
for fn in ./data/klebs/pathogenwatch/assemblies-to-test/*.fasta; do
    sample=$(basename $fn)
    sample="${sample%.*}"
    echo -e "-------------------------------\n"
    echo -e "running mlst on: $sample - $fn"
    
    mlst \
        --threads 2 \ # num of processing threads
        --blastdb $MLST_DB_DIR/blast/mlst.fa \ # specific BLAST-indexed database file containing all MLST allele sequences
        --datadir $MLST_DB_DIR/pubmlst \ # directory containing schem definitions
        --scheme klebsiella \ # specify scheme: gapA, infB, mdh, pgi, phoE, rpoB, and tonB for kleb
        --minid 95 \ # min match identity
        --mincov 10 \ # min match coverage
        --minscore 50 \ # min align score
        ./results/ont/klebsiella/prokka/SRR28370682.fna \ # input
        > ./results/ont/klebsiella/mlst/${sample}.tsv # output
done
```
</details>

##### MLST Output Format

The standard output of MLST analysis is a tabular plain text file with the following columns.

Column | Description
--- | ---
FILE | Input sequence name
SCHEME | The specific bacterial species or genus
ST | Sequence type number
Allelic profile | Depends on how many gene are used in the scheme

> **Note:** When reporting MLST results the scheme used for the profiling must be provided for accurate results interpretation.

<details>
  <summary>
    Click to toggle <b style="color:blue">BIGSdb platform for assigning STs</b>
  </summary>
  <p>
    <a href="https://bigsdb.pasteur.fr/klebsiella" target="_blank">BIGSdb</a> is curated by the Institut Pasteur.
  </p>
</details>

##### MLST Result Interpretation
The most important information of the results in the ST:
- **Know STs:** which matches a database hit a number is assigned.
- **Novel allele combinations:** may be represented as `?` or `novel` if no known database profile is matched
- **Incomplete matches:** may be shown as `ST-like` or with an `*` (asterisk) if most but not all allele profiles match a known database ST
- **Coverage, identity and depth:** are importance to consider in interpratation
- Look out for **STs with specific characteristics**, eg., _K. pneumoniae_ ST258 is a major global clone carrying KPC carbapenemases
- **Clonal complexes:** STs sharing alleles at most loci, usually sharing identical alleles at 5 or more loci; reported as `CC` followed by the number of the central/founding ST
- MLST reuslts can indicate evolutionary relationships

##### Visualising MLST Results
A number of visualisation tools are available, examples:
- [**GrapeTree**](https://achtman-lab.github.io/GrapeTree/MSTree_holder.html): creates hierarchical clustering of MLST data
- [**goeBURST**](https://www.phyloviz.net/goeburst/#Software): a classic tool for visualizing MLST data that focuses on identifying clonal complexes

##### MLST Interpretation Limitations
When analyzing MLST results, be aware of ceratain limitations:

- Limited Resolution: MLST may not distinguish between closely related isolates
- Temporal Dynamics: Does not capture all evolutionary changes over time
- Geographic Bias: Some STs may be overrepresented in databases due to sampling bias
- Horizontal Gene Transfer: May complicate interpretation of evolutionary relationships

### Step 6: Variant Calling and Consensus Assemblies

We will now focus on using an alignment-based comparison approach to identify relationships within pathogen isolates. We begin by retrieving a relevant reference genome. 
```
# Create a ref directory and retrieve reference
mkdir -p ./data/klebs/reference

# Download, unzip ref genome
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/016/305/GCF_000016305.1_ASM1630v1/GCF_000016305.1_ASM1630v1_genomic.gbff.gz \
-O - | gunzip -c > ./data/klebs/reference/GCF_000016305.1_ASM1630v1_genomic.gbff
```

#### Fast Bacterial Variant Calling with Contigs
```
# Create a working directory
mkdir -p ./data/klebs/pathogenwatch/assemblies-to-test
```

```
# Copy assembled genome to the working directory
rsync -avP \
    ./results/ont/klebsiella/dragonflye/SRR28370682.fa \
    ./data/klebs/pathogenwatch/assemblies-to-test/SRR28370682.fasta
```

To model an epidemilogical situation involving different pathogens in circulation, 11 isolate assemblies collected in Kenya between January 14 and January 31, 2019 were retrieved from [PathogenWatch](https://pathogen.watch/genomes/all?country=ke&genusId=570&maxDate=2019-01-31T20%3A59%3A59.999Z&minDate=2018-12-31T21%3A00%3A00.000Z&sort=date&speciesId=573)  as context to ourassembled isolate. They are located in `./data/klebs/pathogenwatch/genomes`. We copy the isolates into the working directory we just created.

```
rsync -avP \
    ./data/klebs/pathogenwatch/genomes/SAMN25722[2,3]{68,97,35,64,03,77}.fasta \
    ./data/klebs/pathogenwatch/assemblies-to-test/
```
Define the temporary directory for `snippy` analysis.

```
export TMPDIR=$(pwd)/results/ont/klebsiella/tmp/snippy/
```
Contig-based variant calling with `snippy`:
```
# Iterate over assemblies
for fn in ./data/klebs/pathogenwatch/assemblies-to-test/*.fasta; do
    sample=$(basename $fn) # retrieve sample name
    sample="${sample%.*}" # remove file extension suffix
    echo -e "-------------------------------\n"
    echo -e "running snippy on: $sample - $fn"
    
    snippy \
        --force \ # overwrite existing dirs
        --prefix $sample \ # output prefix
        --cpus 4 \ # number of processing threads
        --ram 4 \ # memory allocation
        --tmpdir $TMPDIR \ # temporary directory
        --outdir ./results/ont/klebsiella/snippy/$sample \ # output
        --reference ./data/klebs/reference/GCF_000016305.1_ASM1630v1_genomic.gbff \ # reference genome
        --ctgs $fn # assembled contigs
        --report # summarises variants in a human-readable way
    echo -e "-------------------------------\n"
done
```
##### `Snippy` Outputs
Snippy can perform variant calling based on reads or contigs. We used the latter approach which is associated with the following output files (Some of them are universal to either approches except the amount of infomation may differ).

File | Description
--- | ---
`*.tab` | Lists all variants in a human-readable form; favourable for quick epidemilogical insight
`*.vcf` | Universal variant calling format parsible by intergrating bioinformatics tools
`*.gff` | Provides genomic feature context for each variant
`*.txt` | Lists numbers of different variant types
`*.consensus.fa` | FASTA file of consensus sequence; variants are incoporated into the reference
`*.aligned.fa` | The alignment between the reference and the consensus in FASTA format; useful for phylogenetics
`*.bam` | BAM file for **contigs** alignment agains the reference
`*.html` | HTML summary of variants
`*.txt` | A simple summary of the run; lists numbers of different variant types
##### Visualising the Snippy Variants

1. Go to https://igv.org/app/
2. Load the Genome `ref.fa` in the `Genome` tab
3. Load the alignment and its index in the `Tracks` tab
4. For the `ref.fa` select `NC_009648` as chromosome and type `NC_009648:1-200`
   as the region of interest


<br>
<left><img src="img/igv-screenshot.png" alt="Screenshot of IGV" width="1532"/></left>
<br>


##### Build Core and Whole Genome Aligments from Snippy Output

Snippy-core combines multiple Snippy outputs to:

- Identify the conserved-sequence genome ("core genome") which are regions present in all isolates.
- Extract SNPs from these core regions
- Create a multiple sequence alignment of core SNPs (SNP absence/presence-based avoiding the full compilation of variation such as `INS`, `DELS`, variant types)
- Generate various outputs suitable for downstream phylogenetic analysis

#### Run snippy-core
`snippy-core` is a tool that allows you to compile a core **SNP alignment** from multiple Snippy analyses. It's particularly useful for phylogenetic analysis of closely related bacterial isolates.
```
snippy-core \
    --ref ./data/klebs/reference/GCF_000016305.1_ASM1630v1_genomic.gbff \ # ref sequence
    --prefix ./results/ont/klebsiella/snippy-core/core-snp \ #output prefix
    ./results/ont/klebsiella/snippy/* # snippy dir outputs
```
#### Snippy Core Outputs
File | Description
--- | ---
`*.aln` | Alignment of core SNPs in FASTA format
`*.full.aln` | Alignment of entire core genome (including invariant sites)
`*.vcf` | Multi-sample VCF file with all core SNP sites
`*.tab` | Tab-separated table of core SNPs with annotations
`*.text` | Summary statistics text file
`*.ref.fa` | The reference sequence used
`*.aligned.fa` | Contains one sequence per isolate including the reference
<!-- `*.nwk` | Fast neighbor-joining tree based on core SNP alignment -->

#### Cleanup the Snippy SNP Alignment
The resultant `core-snp.full.aln` is loaded with alphabetical encoding (see [Snippy](https://github.com/tseemann/snippy/tree/master) for details) which may not be compatible with downstream analyses like tree-building or recombination-removal. The numerous encodings may be simplified  as below:
```
snippy-clean_full_aln \
    ./results/ont/klebsiella/snippy-core/core-snp.full.aln > \ # snippy-core output
    ./results/ont/klebsiella/snippy-core/core-snp-clean.full.aln # simplified aligment ouput

```

#### Compute Pairwise SNP Distances
`snp-dists` is a basic and fast command-line tool for computing pairwise SNP distances between sequence assemblies. We can use it to quickly generate distance matrices for phylogenetic and epidemiological analyses.
>**Note:** 
`snp-dists` is a stand-alone tool to be installed separately

We can process with phylogentic-analyses in mind--option `-p`:
```
snp-dists -p \ # a phylo-compartible distance matrix
    ./results/ont/klebsiella/snippy-core/core-snp-clean.full.aln > \ # clean alignment
    ./results/ont/klebsiella/snippy-core/core-snp.distance.phylip # distance matrix output
```
Which can allow you generate a phylogentic tree using, for instance `FastTree` (which you must have loaded/installed):
```
FastTree -dist core-snp.distance.phylip > kleb_pneu_phylogenetic_tree.nwk
```
Or default to the absolute matrix distance:
```
snp-dists \ # a phylo-compartible distance matrix
    ./results/ont/klebsiella/snippy-core/core-snp-clean.full.aln > \ # clean alignment
    ./results/ont/klebsiella/snippy-core/core-snp.distance.tsv # distance matrix output
```
### Step 7: Dealing with Recombination
Sometimes you may want to exclude certain regions of the genome from snippy consideration for practical reasons. This is common for genomes like _M.tuberculosis_ where pesky repetitive 
PE/PPE/PGRS genes cause false positives, or masking phage regions.

The tool provides a conveniency by using the `--mask` option with `BED` file of the regions for exclusion.
Any SNPs in those regions will be excluded. 
 
> **Note:** An exclusion `BED` file for _M.tb_ is provided with Snippy in the 
`etc/Mtb_NC_000962.3_mask.bed` folder. 
It is derived from the `XLSX` file from https://gph.niid.go.jp/tgs-tb/.

What about a scenario where the masking `BED` file is not provided?

Gubbins (Genealogies Unbiased By recomBinations In Nucleotide Sequences) is a toolfor detecting and accounting for homologous recombination in bacterial whole genome alignments. Since homologous recombination can obscure the true relatedness/vertical inheritance, it may be desirable to find and mask such regions before phylogenetic inference.

Gubbins:
- Takes a core genome alignment (e.g., from Snippy-core) as input.
- Detects recombination regions using elevated SNP density or maximum likelihood methods (RAxML/FastTree/PhyML).
- Iteratively masks/"remove" recombination and infers a clonal phylogeny.

Predicting recombination regions:
```
# Detect recombination
run_gubbins.py \
    --threads 2 \ # number of processing threads
    --prefix ./results/ont/klebsiella/gubbins/core-snp \ # output prefix
    --iterations 5 \ # max iterations for refining predictions
    --min-snps 3 \ # min number of snp to consider a recombination block
    --min-window-size 100 \ # smalles window size for scanning recomb
    --max-window-size 10000 \# largest window size for scanning recomb
    --filter-percentage 25.0 \ # mask alignment columns with x missing/ambigous info
    ./results/ont/klebsiella/snippy-core/core-snp-clean.full.aln
```

Masking the predicted recombination regions:

```
# Recomb region masking
mask_gubbins_aln.py \
    --aln ./results/ont/klebsiella/snippy-core/core-snp-clean.full.aln \ # original/sbippy alignment
    --gff ./results/ont/klebsiella/gubbins/core-snp.recombination_predictions.gff \ # .gff from gubbins recombinant region annotation
    --out ./results/ont/klebsiella/gubbins/core-snp.masked.aln # output files
```

We can then convert the Gubbins `.gff` masking feature output into a `BED` format file usable with `snippy` as arguement to the `--mask` option.
```
# Convert Gubbins GFF to BED format
python gubbins_to_bed.py \
./results/ont/klebsiella/gubbins/core-snp.recombination_predictions.gff \ # .gff from gubbins recombinant region annotation > ./results/ont/klebsiella/gubbins/gubbins_recomb.bed
```
<details>
    <summary>Click to toggle a <b style='color:blue'>Challenge</b>
    </summary>
Now, can you use the resultant `BED` file to re-run Snippy all the way to building a phylogenrtic tree. How do the trees compare with(out) masking of recombinat regions?
</details>

#### Phylogenetic Analysis of Gubbins Output
We can use another phylogenetic tool--`iqtree`--to visualise the relationships between the recombination-masked isolates from Gubbins.
```
iqtree \
    -m HKY \ # substituion model
    -bb 1000 \ # bootstrap iterations
    -alrt 1000 \ # likelihood ratio test
    -alninfo \ # log aln ino
    -s ./results/ont/klebsiella/gubbins/core-snp.masked.aln \ # input seq aln
    -nt 2 \ # number of threads
    -redo \ # overwrite existing output
    -pre ./results/ont/klebsiella/iqtree/core-snp # prefix for ouputs
```
> **Note:** It is important to use phylogenetic algorithms that take into account SNP alignments. These algorithms usually include some form of ascertainment bias correction that corrects for the 'missing' nucleotides in the alignment that were masked/removed because they did not show polymorphism.

