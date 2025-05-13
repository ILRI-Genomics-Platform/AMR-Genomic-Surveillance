## Login to hpc and create project directory
```
interactive -w compute06 -c 8

cd /var/scratch/$USER

BASEDIR=$(pwd)

echo $BASEDIR


mkdir -p \
$BASEDIR/trainings/ACDC_AMR2025/results/illumina/ecoli/{fastqc,fastp,fastq-scan,shovill,prokka,amrfinder,mlst}


cd $BASEDIR/trainings/ACDC_AMR2025


ln -s /var/scratch/global/jjuma/ACDC_AMR2025/[dpsr]* .
```

## Loading modules
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

## Quality control (QC)
```
fastqc \
    -o ./results/illumina/ecoli/fastqc \
    --noextract \
    -t 2 \
    --nogroup \
    ./data/ecoli/illumina/SRR25008769_1.fastq.gz \
    ./data/ecoli/illumina/SRR25008769_2.fastq.gz


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
    --qualified_quality_phred 25 \
    --unqualified_percent_limit 40 \
    --length_required 20 \
    2> ./results/illumina/ecoli/fastp/SRR25008769.fastp.log

```
## QC Summary Statistics
fastq-scan reads a FASTQ from STDIN and outputs summary statistics (read lengths, per-read qualities, per-base qualities) in JSON format.
```
gzip -cd \
./data/ecoli/illumina/SRR25008769_1.fastq.gz | fastq-scan -g 5249449 > \
./results/illumina/ecoli/fastq-scan/SRR25008769_1-original.json

gzip -cd \
./data/ecoli/illumina/SRR25008769_2.fastq.gz | fastq-scan -g 5249449 > \
./results/illumina/ecoli/fastq-scan/SRR25008769_2-original.json


gzip -cd \
./results/illumina/ecoli/fastp/SRR25008769_R1.trim.fastq.gz | fastq-scan -g 5249449 > \
./results/illumina/ecoli/fastq-scan/SRR25008769_1-trimmed.json

gzip -cd \
./results/illumina/ecoli/fastp/SRR25008769_R2.trim.fastq.gz | fastq-scan -g 5249449 > \
./results/illumina/ecoli/fastq-scan/SRR25008769_2-trimmed.json

```

## De novo assembly pipeline for Illumina paired reads

```
mkdir ./results/illumina/ecoli/tmp_shovill
TMPDIR=./results/illumina/ecoli/tmp_shovill

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
  --cpus 2 \
  --ram 7 \
  --tmpdir $TMPDIR

mv ./results/illumina/ecoli/shovill/SRR25008769/contigs.fa ./results/illumina/ecoli/shovill/SRR25008769/SRR25008769.fa

rm ./results/illumina/ecoli/tmp_shovill
```

## Annotation
```
mkdir ./results/illumina/ecoli/tmp_prokka/
TMPDIR=./results/illumina/ecoli/tmp_prokka/

prokka \
    --evalue 1e-09 \
    --coverage 80 \
    --centre ILRI \
    --cpus 2 \
    --prefix SRR25008769 \
    --locustag SRR25008769 \
    --proteins ./databases/prokka/proteins.faa \
    --force \
    --outdir ./results/illumina/ecoli/prokka \
    ./results/illumina/ecoli/shovill/contigs.fa

rm -r ./results/illumina/ecoli/tmp_prokka/

```


## Full AMRFinderPlus search combining results
```
mkdir ./results/illumina/ecoli/tmp_amrfinder/
TMPDIR=./results/illumina/ecoli/tmp_amrfinder/

AMRFINDER_DB=$(find ./databases/amrfinderplus/2023-11-15.1 -name "AMR.LIB" | sed 's=AMR.LIB==')

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


rm -r ./results/illumina/ecoli/tmp_amrfinder
```

## Multilocus sequence typing
```
mkdir ./databases/mlst/database
tar -xzf ./databases/mlst/mlst.tar.gz -C ./databases/mlst/database
MLST_DB=$(find ./databases/mlst/database/ -name "mlst.fa" | sed 's=blast/mlst.fa==')

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
