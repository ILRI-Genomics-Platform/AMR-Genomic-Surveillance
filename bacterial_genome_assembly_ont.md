



mkdir -p /Users/jjuma/trainings/ACDC_AMR2025/results/ont/klebsiella/{porechop,nanoq,fastq-scan,nanoplot}

# Remove Adapters
porechop \
    --input ./data/klebs/ont/SRR28370682.fastq.gz \
    --format fastq.gz \
    --threads 2 \
    --no_split \
    --output ./results/ont/klebsiella/porechop/SRR28370682_adapter_ont.fastq.gz

# Quality filter
nanoq \
    --min-len 1000 \
    --min-qual 0 \
    --input ./results/ont/klebsiella/porechop/SRR28370682_adapter_ont.fastq.gz \
    --output-type g \
    --output ./results/ont/klebsiella/nanoq/SRR28370682_filt.fastq.gz


# Quality steps before and after QC
gzip -cd \
./data/klebs/ont/SRR28370682.fastq.gz | fastq-scan -g 5249449 > \
./results/ont/klebsiella/fastq-scan/SRR28370682-original.json

gzip -cd \
./results/ont/klebsiella/nanoq/SRR28370682_filt.fastq.gz | fastq-scan -g 5249449 > \
./results/ont/klebsiella/fastq-scan/SRR28370682-final.json


# Nanoplot
NanoPlot \
    --threads 2 \
    --fastq ./data/klebs/ont/SRR28370682.fastq.gz \
    --outdir ./results/ont/klebsiella/nanoplot/ \
    --prefix SRR25570560-original_


<!-- cp results/summary/SRR28370682-original/SRR28370682-original_NanoPlot-report.html
cp results/summary/SRR28370682-original_NanoPlot-report.html -->

NanoPlot \
    --threads 2 \
    --fastq ./results/ont/klebsiella/nanoq/SRR28370682_filt.fastq.gz \
    --outdir ./results/ont/klebsiella/nanoplot/ \
    --prefix SRR28370682-final_

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




