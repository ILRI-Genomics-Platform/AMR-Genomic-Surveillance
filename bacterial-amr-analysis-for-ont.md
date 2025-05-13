
# *K. pneumoniae* AMR analysis workflow using ONT data
---  

###### **_Trainers_**: [John Juma](https://github.com/ajodeh-juma), [Kennedy Mwangi](https://github.com/wanjauk), [Ouso Daniel](https://github.com/ousodaniel) & [Gilbert Kibet](https://github.com/kibet-gilbert)

---

# Project: PRJNA1087001
Benchmark various nanopore-based variant callers on 14 different species.
Samples are sequenced on the latest (September 2023) R10.4.1 Nanopore flowcells and Illumina. Ground truth assemblies are generated for each sample.

https://elifesciences.org/reviewed-preprints/98300

# Set up directories
Before starting the analysis, ensure that you are logged into the HPC, create an interactive session on the assigned compute node, and change directory to the project folder which is `ACDC_AMR2025`.

```
mkdir -p \
results/ont/klebsiella/{porechop,nanoq,fastq-scan,nanoplot,dragonflye,prokka,resfinder,amrfinder,mlst,tmp/{dragonflye,prokka,resfinder,amrfinder,snippy},snippy,snippy-core,gubbins,iqtree}

ln -sf /var/scratch/global/jjuma/ACDC_AMR2025/[dpsr]* .
```


## Retrieve the reference genome in GenBank format

```
mkdir -p ./data/klebs/reference
```

```
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/016/305/GCF_000016305.1_ASM1630v1/GCF_000016305.1_ASM1630v1_genomic.gbff.gz -P  ./data/klebs/reference

```

```
gzip -c -d ./data/klebs/reference/GCF_000016305.1_ASM1630v1_genomic.gbff.gz > ./data/klebs/reference/GCF_000016305.1_ASM1630v1_genomic.gbff
```


# Remove Adapters

```
module load porechop/0.2.4
```

```
porechop \
    --input ./data/klebs/subsampled/SRR28370682.fastq.gz \
    --format fastq.gz \
    --threads 4 \
    --no_split \
    --output ./results/ont/klebsiella/porechop/SRR28370682_adapter_ont.fastq.gz
```

```
module unload porechop/0.2.4
```


# Load modules

```
module load any2fasta/0.4.2
module load nanoq/0.10.0
module load nanoplot/1.42.0
module load bbmap/38.95


module load bwa/0.7.17
module load samtools/1.9
module load racon/1.5.0
module load prodigal/2.6.3  
module load fastp/0.22.0
module load minimap2/2.13
module load spades/3.13.0
module load mlst/2.23.0
module load infernal/1.1.2  
module load fastqc/0.11.9

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

# Quality filter

```
nanoq \
    --min-len 1000 \
    --min-qual 0 \
    --input ./results/ont/klebsiella/porechop/SRR28370682_adapter_ont.fastq.gz \
    --output-type g \
    --output ./results/ont/klebsiella/nanoq/SRR28370682_filt.fastq.gz
```

# Nanoplot

```
NanoPlot \
    --threads 2 \
    --fastq ./data/klebs/subsampled/SRR28370682.fastq.gz \
    --outdir ./results/ont/klebsiella/nanoplot/ \
    --prefix SRR28370682-original_
```

### Copy the report html to local computer from hpc

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


### Copy the report html to local computer from hpc

```
rsync -avP \
    --partial \
    ./results/ont/klebsiella/nanoplot/results/summary/SRR28370682-final_NanoPlot-report.html \
    ~/
```

# *De novo* assembly for ONT reads

### Assemble bacterial isolate genomes from Nanopore reads

### Main steps

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

```
export TMPDIR=./results/ont/klebsiella/tmp/dragonflye/
```

```
dragonflye \
    --reads ./results/ont/klebsiella/nanoq/SRR28370682_filt.fastq.gz \
    --gsize 5248520 \
    --prefix SRR28370682 \
    --outdir ./results/ont/klebsiella/dragonflye \
    --assembler flye \
    --tmpdir ./results/ont/klebsiella/tmp/dragonflye/ \
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


# Assembly evaluation

```
stats.sh in=./results/ont/klebsiella/dragonflye/SRR28370682.fa
```

>**Note**
Unfortunately, the N50 and L50 values generated by stats.sh are switched. N50 should be a length and L50 should be a count. The results table below shows the corrected values based on stats.sh outputs.


# Annotation

```
export TMPDIR=./results/ont/klebsiella/tmp/prokka/
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
```


### Cleanup intermediate files

```
rm -r ./results/ont/klebsiella/prokka/*.pdb ./results/ont/klebsiella/prokka/*.pjs ./results/ont/klebsiella/prokka/*.pot ./results/ont/klebsiella/prokka/*.ptf ./results/ont/klebsiella/prokka/*.pto
```


## Fast bacterial variant calling from NGS reads or contigs

# Convert the GenBank format to Fasta format
```
any2fasta \
    ./data/klebs/reference/GCF_000016305.1_ASM1630v1_genomic.gbff > \
    ./data/klebs/reference/Reference.fasta
```


Here we will use 11 Klebs isolates collected in Kenya between January 14 and January 31, 2019
https://pathogen.watch/genomes/all?country=ke&genusId=570&maxDate=2019-01-31T20%3A59%3A59.999Z&minDate=2018-12-31T21%3A00%3A00.000Z&sort=date&speciesId=573


```
mkdir -p ./data/klebs/pathogenwatch/assemblies-to-test
```

```
rsync -avP --partial \
    ./results/ont/klebsiella/dragonflye/SRR28370682.fa \
    ./data/klebs/pathogenwatch/assemblies-to-test/SRR28370682.fasta
```

```
rsync -avP --partial \
    ./data/klebs/reference/Reference.fasta \
    ./data/klebs/pathogenwatch/assemblies-to-test/Reference.fasta
```

## Select assemblies to test

 <!-- ./data/klebs/pathogenwatch/genomes/SAMN25722[2,3]{68,97,35,64,03,77}.fasta -->

```
rsync -avP --partial \
    /data/klebs/pathogenwatch/genomes/*.fasta
    ./data/klebs/pathogenwatch/assemblies-to-test/
```


# AMR genes detection using ResFinder

```
python -m resfinder \
    -ifa ./results/ont/klebsiella/dragonflye/SRR28370682.fa \
    -o ./results/ont/klebsiella/resfinder/ \
    -s klebsiella \
    --min_cov 0.6 \
    --threshold 0.9 \
    --min_cov_point 0.6 \
    --threshold_point 0.9 \
    --ignore_stop_codons \
    --ignore_indels \
    --acquired \
    --point
```

# Identify AMR and virulence genes in proteins and/or contigs
### Full AMRFinderPlus search combining results

# Full AMRFinderPlus search combining results

Identify acquired antimicrobial resistance genes in bacterial protein and/or assembled nucleotide sequences as well as known resistance-associated point mutations for several taxa. With AMRFinderPlus we added select members of additional classes of genes such as virulence factors, biocide, heat, acid, and metal resistance genes.

>**Note**
AMRFinderPlus reports gene and point mutation presence/absence; it does not
infer phenotypic resistance. Many of the resistance genes detected by
AMRFinderPlus may not be relevant for clinical management or antimicrobial
surveillance. **The AMR genes must be expressed to confer resistance**

AMRFinder can be run in multiple modes with protein sequence as input and/or with DNA sequence as input. To get maximum information it should be run with both protein and nucleotide. When run with protein sequence it uses both BLASTP and HMMER to search protein sequences for AMR genes along with a hierarchical tree of gene families to classify and name novel sequences. With nucleotide sequences it uses BLASTX translated searches and the hierarchical tree of gene families. Adding the `--organism` option enables screening for point mutations in select organisms and suppresses the reporting of some that are extremely common in those organisms.

```
AMRFINDER_DB=$(find /export/apps/amrfinder/4.0.22/data/2025-03-25.1 -name "AMR.LIB" | sed 's=AMR.LIB==')
```
```
export TMPDIR=./results/ont/klebsiella/tmp/amrfinder/
```
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

# AMR detection with ResFinder

```
python -m resfinder \
    -ifa ./results/ont/klebsiella/dragonflye/SRR28370682.fa \
    -o ./results/ont/klebsiella/resfinder/SRR28370682 \
    -s klebsiella \
    --min_cov 0.6 \
    --threshold 0.9 \
    --min_cov_point 0.6 \
    --threshold_point 0.9 \
    --ignore_stop_codons \
    --ignore_indels \
    --acquired \
    --point
```


# Batch AMR detection 

```
for fn in ./data/klebs/pathogenwatch/assemblies-to-test/*.fasta; do
    sample=$(basename $fn)
    sample="${sample%.*}"
    echo -e "-------------------------------\n"
    echo -e "running ResFinder on: $sample - $fn"

    python -m resfinder \
        -ifa $fn \
        -o ./results/ont/klebsiella/resfinder/${sample} \
        -s klebsiella \
        --min_cov 0.6 \
        --threshold 0.9 \
        --min_cov_point 0.6 \
        --threshold_point 0.9 \
        --ignore_stop_codons \
        --ignore_indels \
        --acquired \
        --point
done
```


# MLST 

```
MLST_DB=$(find ./databases/mlst/database/ -name "mlst.fa" | sed 's=blast/mlst.fa==')
```

### Automatic MLST calling from assembled contigs

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

# Assign sequence types

BIGSdb platform curated by the Institute Pasteur
(https://bigsdb.pasteur.fr/klebsiella)


# Batch MLST typing

```
for fn in ./data/klebs/pathogenwatch/assemblies-to-test/*.fasta; do
    sample=$(basename $fn)
    sample="${sample%.*}"
    echo -e "-------------------------------\n"
    echo -e "running mlst on: $sample - $fn"

    mlst \
        --threads 2 \
        --blastdb $MLST_DB/blast/mlst.fa \
        --datadir $MLST_DB/pubmlst \
        --scheme klebsiella \
        --minid 95 \
        --mincov 10 \
        --minscore 50 \
        $fn \
        > ./results/ont/klebsiella/mlst/${sample}.tsv
done
```

# concatenate output

```
cat \
    ./results/ont/klebsiella/mlst/*.tsv > \
    ./results/ont/klebsiella/mlst/klebs-mlst.txt
```

### Running on hpc
```
module purge
module load snippy/4.6.0
```

### Running on local computer in a conda environment

```
conda activate $HOME/miniforge3/envs/snippy-env
```

## SNP calling and phylogenetic analysis

Snippy finds SNPs between a haploid reference genome and your NGS sequence
reads. It will find both substitutions (snps) and insertions/deletions (indels).
It will use as many CPUs as you can give it on a single computer (tested to 64
cores). It is designed with speed in mind, and produces a consistent set of
output files in a single folder. It can then take a set of Snippy results using
the same reference and generate a core SNP alignment (and ultimately a
phylogenomic tree).



```
export TMPDIR=$(pwd)/results/ont/klebsiella/tmp/snippy/
```

```
for fn in ./data/klebs/pathogenwatch/assemblies-to-test/*.fasta; do
    sample=$(basename $fn)
    sample="${sample%.*}"
    echo -e "-------------------------------\n"
    echo -e "running snippy on: $sample - $fn"
    
    snippy \
        --force \
        --prefix $sample \
        --mapqual 60 \
        --basequal 13 \
        --cpus 4 \
        --ram 4 \
        --tmpdir $TMPDIR \
        --outdir ./results/ont/klebsiella/snippy/$sample \
        --reference ./data/klebs/reference/GCF_000016305.1_ASM1630v1_genomic.gbff \
        --ctgs $fn
    echo -e "-------------------------------\n"
done
```

`--mapqual` is the minimum mapping quality to accept in variant calling. 
BWA MEM using `60` to mean a read is "uniquely mapped".

`--basequal` is minimum quality a nucleotide needs to be used in variant
calling. We use `13` which corresponds to error probability of ~5%. It is a
traditional SAMtools value.

The variant calling is done by Freebayes. The key parameters under user control are:

`--mincov` - the minimum number of reads covering a site to be considered (default=10)
`--minfrac` - the minimum proportion of those reads which must differ from the reference
`--minqual` - the minimum VCF variant call "quality" (default=100)


>**Note**
If you run Snippy with the `--report` option it will automatically run `snippy-vcf_report` and generate a `snps.report.txt`

# Visualize the variants

1. Go to https://igv.org/app/
2. Load the Genome `ref.fa` in the `Genome` tab
3. Load the alignment and its index in the `Tracks` tab
4. For the `ref.fa` select `NC_009648` as chromosome and type `NC_009648:1-200`
   as the region of interest


<br>
<left><img src="img/igv-screenshot.png" alt="Screenshot of IGV" width="1532"/></left>
<br>


# Build whole+core genome aligment from Snippy folders

If you call SNPs for multiple isolates from the same reference, you can produce 
an alignment of "core SNPs" which can be used to build a high-resolution 
phylogeny (ignoring possible recombination). A "core site" is a genomic position Q
that is present in all the samples. A core site can have the same nucleotide in 
every sample ("monomorphic") or some samples can be different 
("polymorphic" or "variant"). If we ignore the complications of "ins", "del" 
variant types, and just use variant sites, these are the "core SNP genome".

If you want to mask certain regions of the genome, you can provide a BED file 
with the `--mask` parameter. Any SNPs in those regions will be excluded. 
This is common for genomes like M.tuberculosis where pesky repetitive 
PE/PPE/PGRS genes cause false positives, or masking phage regions. 
A `--mask` bed file for M.tb is provided with Snippy in the 
`etc/Mtb_NC_000962.3_mask.bed` folder. 
It is derived from the XLSX file from https://gph.niid.go.jp/tgs-tb/

# Run snippy-core
```
snippy-core \
    --maxhap 100 \
    --mask-char N \
    --ref ./data/klebs/reference/GCF_000016305.1_ASM1630v1_genomic.gbff \
    --prefix ./results/ont/klebsiella/snippy-core/core-snp \
    ./results/ont/klebsiella/snippy/*
```

## Output
`.aln` - A core SNP alignment in the `--aformat` format (default FASTA)
`.full.aln` - A whole genome SNP alignment (includes invariant sites)


## Cleanup the alignment
You can remove all the "weird" characters and replace them with N using the 
included snippy-clean_full_aln. This is useful when you need to pass it to a tree-building or recombination-removal tool:

```
snippy-clean_full_aln \
    ./results/ont/klebsiella/snippy-core/core-snp.full.aln > \
    ./results/ont/klebsiella/snippy-core/core-snp-clean.full.aln
```

## Compute snp distances

>**Note** 
Requires `snp-dists`

```
snp-dists \
    ./results/ont/klebsiella/snippy-core/core-snp-clean.full.aln > \
    ./results/ont/klebsiella/snippy-core/core-snp.distance.tsv
```


# Detect recombination

```
run_gubbins.py \
    --threads 8 \
    --prefix ./results/ont/klebsiella/gubbins/core-snp \
    --iterations 5 \
    --min-snps 3 \
    --min-window-size 100 \
    --max-window-size 10000 \
    --filter-percentage 25.0 \
    ./results/ont/klebsiella/snippy-core/core-snp-clean.full.aln
```

## Create masked alignment

```
mask_gubbins_aln.py \
    --aln ./results/ont/klebsiella/snippy-core/core-snp-clean.full.aln \
    --gff ./results/ont/klebsiella/gubbins/core-snp.recombination_predictions.gff \
    --out ./results/ont/klebsiella/gubbins/core-snp.masked.aln
```


# Phylogenetic analysis

>**Note**
Use phylogenetic algorithms that take into account SNP alignments. These
algorithms usually include some form of ascertainment bias correction that
corrects for the 'missing' nucleotides in the alignment that were removed
because they did not show polymorphism.

If working with a recombining species, the alignment will also include SNPs
introduced through recombination, unless recombination detection and masking was
previously performed.

```
iqtree \
    -m HKY \
    -bb 1000 \
    -alrt 1000 \
    -wbt \
    -wbtl \
    -alninfo \
    -s ./results/ont/klebsiella/gubbins/core-snp.masked.aln \
    -nt 8 \
    -ntmax 8 \
    -redo \
    -pre ./results/ont/klebsiella/iqtree/core-snp
```