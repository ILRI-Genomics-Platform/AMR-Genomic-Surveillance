# Set up working Dir:
```
mkdir -p /var/scratch/$USER  
cd /var/scratch/$USER  
```

```
mkdir -p AMR-Genomic-Surveillance/{data,scripts}
cd AMR-Genomic-Surveillance/  
```
# Data  
```
mkdir -p ./data/ont/{klebs,ecoli}
ln -s /var/scratch/global/jjuma/ACDC_AMR2025/data/klebs/ont/* ./data/ont/klebs/
ln -s /var/scratch/global/jjuma/ACDC_AMR2025/data/ecoli/ont/* ./data/ont/ecoli/
```
## Alternative 2: Download data from NCBI/ENA:  
SRA: https://www.ncbi.nlm.nih.gov/sra; search: `klebsiella pneumoniae`
checkboxes to check in SRA:  
`Type`: `genome`; `Platform`: `Oxford Nanopore` OR `Illumina`; `File Type`: `fastq`  
Note: With SRA, first download SRA general format using `wget`, then convert SRA format to FASTQ using `fastq-dump`  

Alternatively download the `fastq.gz` format directly from  the European Nucleotide Archive (ENA) - which mirrors many datasets in NCBI's SRA. Search the accession identified from SRA above in ENA site: https://www.ebi.ac.uk/ena/browser/home  

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/024/SRR21465924/SRR21465924_1.fastq.gz -P ./data/ont/klebs/
```
# on hpc - load modules  

```
module load bbmap/38.95  
module load bwa/0.7.17  
module load samtools/1.9  
module load racon/1.5.0  
module load unicycler/0.4.7  
module load prodigal/2.6.3  
module load fastp/0.22.0  
module load porechop/0.2.4  
module load minimap2/2.13  
module load spades/3.13.0  
module load mlst/2.23.0  
module load infernal/1.1.2  
module load fastqc/0.11.9  
module load any2fasta/0.4.2  
module load medaka/0.8.2  
module load velvet/1.2.10  
module load hmmer/3.3  
module load prokka/1.14.6   
module load lighter/1.1.2  
module load flye/2.4.2  
module load megahit/1.2.9  
module load bowtie2/2.3.4.1  
module load bedtools/2.29.0  
module load flash/1.2.11  
module load htslib/1.9  
module load miniasm/0.3  
module load blast/2.7.1+  
module load barrnap/0.9  
```
Create analysis directories:
```
mkdir -p ~/trainings/ACDC_AMR2025/results/ont/klebsiella/{porechop,nanoq,fastq-scan,nanoplot,dragonflye,prokka,amrfinder,mlst}
```

# Remove Adapters
```
porechop \
    --input ./data/klebs/ont/SRR28370682.fastq.gz \
    --format fastq.gz \
    --threads 2 \
    --no_split \
    --output ./results/ont/klebsiella/porechop/SRR28370682_adapter_ont.fastq.gz
```

# Quality filter
```
nanoq \
    --min-len 1000 \
    --min-qual 0 \
    --input ./results/ont/klebsiella/porechop/SRR28370682_adapter_ont.fastq.gz \
    --output-type g \
    --output ./results/ont/klebsiella/nanoq/SRR28370682_filt.fastq.gz
```

fastq-scan reads a FASTQ from STDIN and outputs summary statistics (read lengths, per-read qualities, per-base qualities) in JSON format.

# Quality steps before and after QC
```
gzip -cd \
./data/klebs/ont/SRR28370682.fastq.gz | fastq-scan -g 5248520 > \
./results/ont/klebsiella/fastq-scan/SRR28370682-original.json
```
```
gzip -cd \
./results/ont/klebsiella/nanoq/SRR28370682_filt.fastq.gz | fastq-scan -g 5248520 > \
./results/ont/klebsiella/fastq-scan/SRR28370682-final.json
```

# Nanoplot
```
NanoPlot \
    --threads 2 \
    --fastq ./data/klebs/ont/SRR28370682.fastq.gz \
    --outdir ./results/ont/klebsiella/nanoplot/ \
    --prefix SRR28370682-original_
```
# copy the report html to local computer from hpc
```
rsync -avP \
    --partial \
    ./results/ont/klebsiella/nanoplot/results/summary/SRR28370682-original_NanoPlot-report.html \
    ~/
```
```
NanoPlot \
    --threads 2 \
    --fastq ./results/ont/klebsiella/nanoq/SRR28370682_filt.fastq.gz \
    --outdir ./results/ont/klebsiella/nanoplot/ \
    --prefix SRR28370682-final_
```
# Tools
```
- bbduk: 39.19
- fastp: 0.24.0
- fastqc: 0.12.1
- fastq-scan: 1.0.1
- lighter: 1.1.3
- nanoplot: 1.44.1
- nanoq: 0.10.0
- pigz: 2.8
- porechop: 0.2.4
- rasusa: 2.1.0
```

# assembly
```
dragonflye \
    --reads ./results/ont/klebsiella/nanoq/SRR28370682_filt.fastq.gz \
    --gsize 5248520 \
    --prefix SRR28370682 \
    --outdir ./results/ont/klebsiella/dragonflye \
    --assembler flye \
    --tmpdir ./results/ont/klebsiella/ \
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
```

# Tools
```
- any2fasta: 0.4.2
- assembly-scan: 1.0.0
- bwa: 0.7.19-r1273
- dragonflye: 1.2.1
- flash: 1.2.11
- flye: 2.9.5-b1801
- medaka: 2.0.1
- megahit: 1.2.9
- miniasm: 0.3-r179
- minimap2: 2.28-r1209
- nanoq: 0.10.0
- pigz: 2.8
- pilon: 1.24
- racon: 1.5.0
- rasusa: 2.1.0
- raven: 1.8.3
- samclip: 0.4.0
- samtools: 1.21
- shovill: 1.1.0
- shovill-se: 1.1.0
- skesa: 2.5.1
- spades.py: 4.1.0
- velvetg: 1.2.10
- velveth: 1.2.10
- unicycler: 0.5.1
```

# annotation
```
mkdir ./results/ont/klebsiella/tmp_prokka/
TMPDIR=./results/ont/klebsiella/tmp_prokka/
```
```
prokka \
    --evalue 1e-09 \
    --coverage 80 \
    --centre ILRI \
    --cpus 2 \
    --prefix SRR28370682 \
    --locustag SRR28370682 \
    --proteins ./databases/prokka/proteins.faa \
    --force \
    --outdir ./results/ont/klebsiella/prokka \
    ./results/ont/klebsiella/dragonflye/SRR28370682.fa

rm -r ./results/ont/klebsiella/tmp_prokka/
```
# Make blastdb of contigs, genes, proteins
```
mkdir ./results/ont/klebsiella/prokka/blastdb
```
```
cat ./results/ont/klebsiella/prokka/SRR28370682.fna | \
    makeblastdb \
    -dbtype "nucl" \
    -title "Assembled contigs for SRR28370682" \
    -out ./results/ont/klebsiella/prokka/blastdb/SRR28370682.fna
```
```
cat ./results/ont/klebsiella/prokka/SRR28370682.ffn | \ 
    makeblastdb \
    -dbtype "nucl" \
    -title "Predicted genes sequences for SRR28370682" \
    -out ./results/ont/klebsiella/prokka/blastdb/SRR28370682.ffn
```
```
cat ./results/ont/klebsiella/prokka/SRR28370682.faa | \
    makeblastdb \
    -dbtype "prot" \
    -title "Predicted protein sequences for SRR28370682" \
    -out ./results/ont/klebsiella/prokka/blastdb/SRR28370682.faa
```
# Compress
```
tar -cvf - blastdb/ | gzip -c > SRR28370682/SRR28370682-blastdb.tar.gz
```
```
gzip SRR28370682/*.gff
gzip SRR28370682/*.gbk
gzip SRR28370682/*.fna
gzip SRR28370682/*.faa
gzip SRR28370682/*.ffn
gzip SRR28370682/*.sqn
gzip SRR28370682/*.fsa
gzip SRR28370682/*.tbl -->
```
# Cleanup intermediate files
```
rm -r ./results/ont/klebsiella/prokka/*.pdb ./results/ont/klebsiella/prokka/*.pjs ./results/ont/klebsiella/prokka/*.pot ./results/ont/klebsiella/prokka/*.ptf ./results/ont/klebsiella/prokka/*.pto
```

# Tools
```
 - makeblastdb: 2.16.0+
 - prokka: 1.14.6
```

# Identify AMR and virulence genes in proteins and/or contigs and print a report
```
AMRFINDER_DB=$(find ./databases/amrfinderplus/2023-11-15.1 -name "AMR.LIB" | sed 's=AMR.LIB==')
```

# Full AMRFinderPlus search combining results
```
amrfinder \
    --nucleotide ./results/ont/klebsiella/prokka/SRR28370682.fna \
    --protein ./results/ont/klebsiella/prokka/SRR28370682.faa \
    --gff ./results/ont/klebsiella/prokka/SRR28370682.gff \
    --annotation_format prokka \
    --organism Klebsiella_pneumoniae \
    --plus \
    --ident_min -1 \
    --coverage_min 0.5 \
    --translation_table 11 \
    --database $AMRFINDER_DB \
    --threads 2 \
    --name SRR28370682 > ./results/ont/klebsiella/amrfinder/SRR28370682.tsv
```
# Tools
```
 - amrfinderplus: 4.0.19
 - amrfinderplus-database: 4.0.19
```

# MLST 
```
mkdir ./databases/mlst/database
tar -xzf ./databases/mlst/mlst.tar.gz -C ./databases/mlst/database
MLST_DB=$(find ./databases/mlst/database/ -name "mlst.fa" | sed 's=blast/mlst.fa==')
```

# Automatic MLST calling from assembled contigs
```
mlst \
    --threads 2 \
    --blastdb $MLST_DB/blast/mlst.fa \
    --datadir $MLST_DB/pubmlst \
    --scheme klebsiella \
    --minid 95 \
    --mincov 10 \
    --minscore 50 \
    ./results/ont/klebsiella/prokka/SRR28370682.fna \
    > ./results/ont/klebsiella/mlst/SRR28370682.tsv
```
# Tools
```
- mlst: 2.23.0
- mlst-database: 2.23.0-20240325
```
