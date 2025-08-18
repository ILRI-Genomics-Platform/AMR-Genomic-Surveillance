
# Bioinformatics Analysis for Antimicrobial Resistance Genomic Surveillance: Illumina data; *Escherichia coli*
---  

###### **_Trainers_**: [John Juma](https://github.com/ajodeh-juma), [Kennedy Mwangi](https://github.com/wanjauk), [Ouso Daniel](https://github.com/ousodaniel) & [Gilbert Kibet](https://github.com/kibet-gilbert)

---

## Introduction  

Escherichia coli (*E. coli*) is a **gram-negative, facultative anaerobic, rod-shaped, coliform bacterium**. It is **mostly harmless** and can be found in various environments such as **soil, water, vegetables, undercooked meats, and milk**.  

In humans, *E. coli* is a **common component of the gut microbiota**, aiding in the **synthesis of vitamins**. However, some strains are **pathogenic**, particularly in **infants, the elderly, and immunocompromised individuals**, leading to **intestinal and extraintestinal diseases** such as:  
- **Urinary tract infections (UTI)**
- **Pneumonia**
- **Bacteremia**
- **Peritonitis**
- **Nosocomial infections**

Additionally, *E. coli* has demonstrated **resistance to beta-lactamase antibiotics** (e.g., **cephalosporins, monobactams**) and **carbapenems** (e.g., **imipenem, ertapenem, meropenem**).

You can read more about it's genome structure, diversity, strains (phylogroups, pathogenicity and serotypes) here: [Introduction to *E. coli*](Intro_Escherichia_coli.md)

---

## Scope of the Tutorial

By now we have done whole-genome assembly of *K. pneumoniae* based on long read ONT data. Now we will explore the alternative of using short read Illumina data. For this we will focus on *E. coli*. We will explore the tools and steps used in such analysis and ultimately conduct AMR analysis.


## Setting up the Bioinformatics analysis environment

### Set up directories  

 - Before starting the analysis, ensure that you are logged into the HPC, create an interactive session on the assigned compute node, and change directory to the project folder which is `ACDC_AMR2025`.  
 - First log in to the HPC.   

```
ssh <user_name>@hpc.ilri.cgiar.org
```
 - Now let us secure a four of CPUs in one of the HPC nodes.  

If your username (`user**`) ends with an ***Odd Number*** (1,3,5,7,9) use `compute05` and if it ends with n ***even number*** (2,4,6,8,0) use `compute06`.   
>Compute05  
```
interactive -w compute05 -c 2 -J amr-surveillance -p batch
```
>Compute06  
```
interactive -w compute06 -c 2 -J amr-surveillance -p batch
```
 - Setting up the project directory structure.  

We will be conducting our a anlysis from a directory in the `scratch` space of the HPC. It is a temporary storage area designed for high-speed data access and short-term file storage. It is cleaned continually (say every night) of files older than *X* time (90 days).

 - Ensure the course directory is created and change into it:

```bash
mkdir -p /var/scratch/$USER/ACDC_AMR2025
cd /var/scratch/$USER/ACDC_AMR2025
```

 - Within our project directory let us create a structured project directories:

```bash
mkdir -p \
results/illumina/ecoli/{fastqc,fastp,fastq-scan,shovill,prokka,resfinder,amrfinder,mlst,tmp/{shovill,prokka,resfinder,amrfinder}}
```

 - Then let us link some data.  
 
Our analysis data is in the path `/var/scratch/global/jjuma/ACDC_AMR2025/`. We will create a symbolic link of the folders with data to current directory.

```bash
ln -sf /var/scratch/global/jjuma/ACDC_AMR2025/[dps]* .
```

### Load modules

```
module load fastqc/0.11.9
module load fastp/0.22.0
module load bbmap/38.95
module load shovill/1.1.0
module load prokka/1.14.6
module load resfinder/4.6.0
module load mlst/2.23.0
module load amrfinder/4.0.22
```

## Bioinformatics Analysis

- If by any chance you wish to explore other `E. coli` datasets from SRA or ENA:

<details>
    <summary>
        Click to toggle contents of 
        <b style='color: blue'> Downloading data from SRA or ENA
        </b>
    </summary>
    You can access public *E. coli* data directly from SRA or ENA as follows:

### Download data from NCBI/ENA:  
### Option 1: SRA
 - Go to SRA (https://www.ncbi.nlm.nih.gov/sra)  
 - search: `Escherichia coli 0157`  
 - checkboxes to check in SRA:  
`   - Type`: `genome`;  
    - `Library Layout`: `paired`;
    - `Platform`: `Illumina`;
    - `File Type`: `fastq` 
    - Based on the description, identify a good dataset e.g there is one that is described as "Sequencing of E. coli resistant to 3rd generation cephalosporins isolated from blood culture"; SRA accession: `SRR32302053`. Click on it.

Note: With SRA, first download SRA general format using `wget`, then convert SRA format to FASTQ using `fastq-dump`

### Option 2: ENA
 - Download the fastq.gz directly from  the European Nucleotide Archive (ENA) - which mirrors many datasets in NCBI's SRA. Search the accession identified from SRA above in ENA site: https://www.ebi.ac.uk/ena/browser/home

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR323/053/SRR32302053/SRR32302053_1.fastq.gz ./
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR323/053/SRR32302053/SRR32302053_2.fastq.gz ./
```
</details>


### Assessing Read Quality using fastqc before quality trimming  

The **raw `FASTQ`** files have sequences as generated by the sequencer. They include **poor quality reads**, **adapters**, **PhiX** and some reads may be **duplicates**. We need to check the quality of these suquences and clean up if we need to.  
We will use our first module - **`FASTQC`**, a bioinformatics software used to analyse quality of raw FASTQ format reads and generate visual plots. The report generated by this step will have information on **Number of reads**,**Sequence Quality**, **Sequence length**, **GC content**, **Adapter content**, **duplication rate** and others. Run the command below to execute fastqc on our reads:

```bash
fastqc \
    -o ./results/illumina/ecoli/fastqc \
    --noextract \
    -t 2 \
    --nogroup \
    ./data/ecoli/illumina/SRR25008769_1.fastq.gz \
    ./data/ecoli/illumina/SRR25008769_2.fastq.gz
```
Takes less than 1 minute

 - Now download the results of the fastqc command to your local laptop for evaluation. The results are in a `HTML` file.  
 - First copy to `~/`(home).

```bash
cp ./results/illumina/ecoli/fastqc/*html ~/
```
 - Now copy the report to a directory on **`local pc`**

```bash
rsync -avP --partial <USERXX>@hpc.ilri.cgiar.org:~/SRR25008769*.html ~/
```
 - Open the HTML file and explore.
 - Links available online: [SRR25008769_1_fastqc.html](https://hpc.ilri.cgiar.org/~gkibet/AMR-Genomic-Surveillance/SRR25008769_1_fastqc.html); [SRR25008769_2_fastqc.html](https://hpc.ilri.cgiar.org/~gkibet/AMR-Genomic-Surveillance/SRR25008769_2_fastqc.html)  

### Quality Trimming fastq files with fastp and Trims adapter sequences

After assessing the quality, we will proceed and do some Quality Control (QC). With **`fastp`** module we will trim reads with qualities less than 20 phred score, remove adapters, remove duplicates and remove PhiX if any.

```bash
fastp \
    --in1 ./data/ecoli/illumina/SRR25008769_1.fastq.gz \
    --in2 ./data/ecoli/illumina/SRR25008769_2.fastq.gz \
    --out1 ./results/illumina/ecoli/fastp/SRR25008769_R1.trim.fastq.gz \
    --out2 ./results/illumina/ecoli/fastp/SRR25008769_R2.trim.fastq.gz \
    --detect_adapter_for_pe \
    --thread 2 \
    --json ./results/illumina/ecoli/fastp/SRR25008769.fastp.json \
    --html ./results/illumina/ecoli/fastp/SRR25008769.fastp.html \
    --cut_mean_quality 20 \
    --cut_front \
    --cut_tail \
    --cut_window_size 4 \
    --qualified_quality_phred 25 \
    --unqualified_percent_limit 40 \
    --length_required 20 \
    2>&1 | tee ./results/illumina/ecoli/fastp/SRR25008769.fastp.log
```
Takes less than 1 minute

 - `fastp` also generates a report in `HTML` format. Let us download it and explore. First copy it to home (`~/`)

```bash
cp ./results/illumina/ecoli/fastp/SRR25008769.fastp.html ~/
```

 - Then download:

```bash
rsync -avP --partial <USERXX>@hpc.ilri.cgiar.org:~/SRR25008769.fastp.html ~/
```
 - Examine the report.
 - Links available online: [SRR25008769.fastp.html](https://hpc.ilri.cgiar.org/~gkibet/AMR-Genomic-Surveillance/SRR25008769.fastp.html)

### De novo assembly pipeline for Illumina paired reads
Assemble bacterial isolate genomes from Illumina paired-end reads

Shovill is a pipeline which uses SPAdes at its core, but alters the steps before
and after the primary assembly step to get similar results in less time. Shovill
also supports other assemblers like SKESA, Velvet and Megahit, so you can take
advantage of the pre- and post-processing the Shovill provides with those too.


>**Note**
Shovill is for isolate data only, primarily small haploid organisms. 
It will NOT work on metagenomes or larger genomes. 
Please use Megahit directly instead.


#### Main steps

1. Estimate genome size and read length from reads (unless --gsize provided) - mash
2. Reduce FASTQ files to a sensible depth (default --depth 100) - seqtk
3. Trim adapters from reads (with --trim only) - trimmomatic
4. Conservatively correct sequencing errors in reads - lighter
5. Pre-overlap ("stitch") paired-end reads - flash
6. Assemble with SPAdes/SKESA/Megahit with modified kmer range and PE + long SE reads
7. Correct minor assembly errors by mapping reads back to contigs (bwa-mem + pilon)
8. Remove contigs that are too short, too low coverage, or pure homopolymers
9. Produce final FASTA with nicer names and parseable annotations


```
export TMPDIR="./results/illumina/ecoli/tmp/shovill"
```

```
shovill \
  --R1 ./results/illumina/ecoli/fastp/SRR25008769_R1.trim.fastq.gz \
  --R2 ./results/illumina/ecoli/fastp/SRR25008769_R2.trim.fastq.gz \
  --gsize 5249449 \
  --outdir ./results/illumina/ecoli/shovill/SRR25008769 \
  --assembler skesa \
  --minlen 500 \
  --mincov 2 \
  --force \
  --keepfiles \
  --depth 0 \
  --noreadcorr \
  --namefmt "SRR25008769_%05d" \
  --cpus 4 \
  --ram 16 \
  --tmpdir $TMPDIR
```

```
mv ./results/illumina/ecoli/shovill/SRR25008769/contigs.fa ./results/illumina/ecoli/shovill/SRR25008769/SRR25008769.fa
```

### Assembly evaluation

```
stats.sh in=./results/illumina/ecoli/shovill/SRR25008769/SRR25008769.fa
```

>**Note**
Unfortunately, the N50 and L50 values generated by stats.sh are switched. N50
should be a length and L50 should be a count. The results table below shows the
corrected values based on stats.sh outputs.


### Annotation
```
export TMPDIR=./results/illumina/ecoli/tmp/prokka/
```

```
prokka \
    --evalue 1e-09 \
    --centre ILRI \
    --coverage 80 \
    --cpus 2 \
    --prefix SRR25008769 \
    --locustag SRR25008769 \
    --proteins ./databases/prokka/proteins.faa \
    --force \
    --outdir ./results/illumina/ecoli/prokka \
    ./results/illumina/ecoli/shovill/SRR25008769/SRR25008769.fa
```

### AMR genes detection using ResFinder

```
python -m resfinder \
    -ifa ./results/illumina/ecoli/shovill/SRR25008769/SRR25008769.fa \
    -o ./results/illumina/ecoli/resfinder/ecoli \
    -s ecoli \
    --min_cov 0.6 \
    --threshold 0.9 \
    --min_cov_point 0.6 \
    --threshold_point 0.9 \
    --ignore_stop_codons \
    --ignore_indels \
    --acquired \
    --point
```

### Full AMRFinderPlus search combining results

Identify acquired antimicrobial resistance genes in bacterial protein and/or assembled nucleotide sequences as well as known resistance-associated point mutations for several taxa. With AMRFinderPlus we added select members of additional classes of genes such as virulence factors, biocide, heat, acid, and metal resistance genes.

>**Note**
AMRFinderPlus reports gene and point mutation presence/absence; it does not
infer phenotypic resistance. Many of the resistance genes detected by
AMRFinderPlus may not be relevant for clinical management or antimicrobial
surveillance. **The AMR genes must be expressed to confer resistance**

AMRFinder can be run in multiple modes with protein sequence as input and/or with DNA sequence as input. To get maximum information it should be run with both protein and nucleotide. When run with protein sequence it uses both BLASTP and HMMER to search protein sequences for AMR genes along with a hierarchical tree of gene families to classify and name novel sequences. With nucleotide sequences it uses BLASTX translated searches and the hierarchical tree of gene families. Adding the `--organism` option enables screening for point mutations in select organisms and suppresses the reporting of some that are extremely common in those organisms.


```
export TMPDIR=./results/illumina/ecoli/tmp/amrfinder/
```

```
AMRFINDER_DB=$(find /export/apps/amrfinder/4.0.22/data/2025-03-25.1 -name "AMR.LIB" | sed 's=AMR.LIB==')
```

Purge all the loaded modules and load `amrfinder` and `mlst`.
```
module purge
module load mlst/2.23.0
module load amrfinder/4.0.22
```

```
amrfinder \
    --nucleotide ./results/illumina/ecoli/prokka/SRR25008769.fna \
    --protein ./results/illumina/ecoli/prokka/SRR25008769.faa \
    --gff ./results/illumina/ecoli/prokka/SRR25008769.gff \
    --annotation_format prokka \
    --organism Escherichia \
    --plus \
    --ident_min -1 \
    --coverage_min 0.5 \
    --translation_table 11 \
    --database $AMRFINDER_DB \
    --threads 2 \
    --name SRR25008769 > ./results/illumina/ecoli/amrfinder/SRR25008769.tsv
```

### Multilocus sequence typing

```
MLST_DB=$(find /export/apps/mlst/2.23.0/db -name "mlst.fa" | sed 's=blast/mlst.fa==')
```

```
mlst \
    --threads 2 \
    --blastdb $MLST_DB/blast/mlst.fa \
    --datadir $MLST_DB/pubmlst \
    --scheme ecoli_achtman_4 \
    --minid 95 \
    --mincov 10 \
    --minscore 50 \
    ./results/illumina/ecoli/prokka/SRR25008769.fna \
    > ./results/illumina/ecoli/mlst/SRR25008769.tsv
```

