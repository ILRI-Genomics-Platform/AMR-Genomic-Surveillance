
# *Monkeypox virus* analysis workflow using Illumina data
---  

###### **_Trainers_**: [John Juma](https://github.com/ajodeh-juma), [Kennedy Mwangi](https://github.com/wanjauk), [Ouso Daniel](https://github.com/ousodaniel) & [Gilbert Kibet](https://github.com/kibet-gilbert)

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
      - [Step 1: Load required modules](#step-1-load-required-modules)
      - [Step 2: Create a virtual environment and install dependencies](#step-2-create-a-virtual-environment-and-install-dependencies)
      - [Step 3: Remove human reads](#step-3-remove-human-reads)
      - [Step 4: Trim adapter sequences](#step-4-trim-adapter-sequences)
      - [Step 5: Get primer scheme](#step-5-get-primer-scheme)
      - [Step 6: Alignment of reads](#step-6-alignment-of-reads)
        - [Index the reference genome](#index-the-reference-genome)
        - [Align reads to the reference genome](#align-reads-to-the-reference-genome)
        - [Sort alignment](#sort-alignment)
        - [Index alignment](#index-alignment)
      - [Step 7: Trim alignments from an amplicon scheme](#step-7-trim-alignments-from-an-amplicon-scheme)
      - [Step 8: Call variants](#step-9-call-variants)
        - [Make depth mask, split variants into ambiguous/consensus](#make-depth-mask-split-variants-into-ambiguous-consensus)
        - [Normalize variant records into canonical VCF representation](#normalize-variant-records-into-canonical-vcf-representation)
      - [Step 9: Generate consensus genome](#step-9-generate-consensus-genome)
        - [Apply ambiguous variants first using IUPAC codes](#apply-ambiguous-variants-first-using-iupac-codes)
        - [Get viral contig name from reference](#get-viral-contig-name-from-reference)
        - [Apply remaning variants, including
          indels](#apply-remaning-variants-including-indels)
        - [Coverage metrics and visualization with IGV](#coverage-metrics-and-visualization-with-igv)
      - [Step 10: Squirrel - Some QUIck Reconstruction to Resolve Evolutionary Links](#step-10-squirrel-some-quick-reconstruction-to-resolve-evolutionary-links)
        - [Add additional sequences retrieved from Pathoplexus](#add-additional-sequences-retrieved-from-pathoplexus)
      

## Set-up  
### Workshop Environmnet
The workshop is mainly based on working in a HPC environment. However, local
set-ups will be provisioned, although some steps of the workflows may be
impossible to execute on laptops. For access to the HPC, you will use your
personal computers to log into the ILRI HPC cluster, which operates on a Linux
operating system. If using a non-Unix operating system you will require a
program that enables interfacing with the Linux HPC environment.


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
mkdir -p /var/scratch/$USER/ACDC_AMR2025
```

# Change directory into the created directory

```
cd /var/scratch/$USER/ACDC_AMR2025
```

>**Note** Once inside the `hpc`, all instances of ```$USER``` will be equivalent to the `hpc` username that you were assigned. Your username, by default, is stored in a variable called `USER`. By using it, you will not have to type your username, rather, your shell will automatically retrieve your username which is the value stored in the `USER` variable. The `$` (dollar) prefix is used to retrieve the value of that variable.

2. *Create project sub-directories:*

```
mkdir -p \
results/mpox/{fastqc,fastp,hostile,trim_galore,primerschemes,bwa/{index,alignment},tmp/squirrel,primertrimmed,primer-trimmed,freebayes,pathoplexus-data,all_consensus,squirrel}
```

3. *Create symblic links to the required resources

```
ln -sf /var/scratch/global/jjuma/ACDC_AMR2025/[dpsr]* . 
```


4. *Retrieve files from the GitHub repository

```
mkdir -p ./utils
```

<!-- ```
mkdir -p ./utils
wget -c https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/artic-mpox-environment.yml -P ./utils
``` -->

<!-- ```
cp /var/scratch/global/jjuma/ACDC_AMR2025/artic-requirements.txt . 
``` -->

```
wget -c https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/artic-requirements.txt -P ./utils
```

### Step 1: Load required modules
```
module load hostile/2.0.0
module load fastp/0.22.0
module load seqtk/1.3
module load bwa/0.7.17
module load freebayes/1.3.4
module unload bcftools/1.17
module load bcftools/1.13
module load squirrel/1.1.2
```

### Step 2: Create a virtual and install dependencies

```
python3 -m venv ./py3env
```

```
source ./py3env/bin/activate
```

```
python3 -m pip install -r ./utils/artic-requirements.txt
```

### Step 3: Remove human reads

```
hostile clean \
    --fastq1 ./data/mpox/illumina/SRR21755837_1.fastq.gz \
    --fastq2 ./data/mpox/illumina/SRR21755837_2.fastq.gz \
    --output ./results/mpox/hostile/ \
    --threads 4 \
    --index ./databases/hostile/human-t2t-hla
```

### Step 4: Trim adapter sequences

```
fastp \
    --in1 ./results/mpox/hostile/SRR21755837_1.clean_1.fastq.gz \
    --in2 ./results/mpox/hostile/SRR21755837_2.clean_2.fastq.gz \
    --out1 ./results/mpox/fastp/SRR21755837_R1.trim.fastq.gz \
    --out2 ./results/mpox/fastp/SRR21755837_R2.trim.fastq.gz \
    --detect_adapter_for_pe \
    --thread 4 \
    --json ./results/mpox/fastp/SRR21755837.fastp.json \
    --html ./results/mpox/fastp/SRR21755837.fastp.html \
    --cut_mean_quality 20 \
    --qualified_quality_phred 25 \
    --unqualified_percent_limit 40 \
    --length_required 20 \
    2> ./results/mpox/fastp/SRR21755837.fastp.log
```
 
### Step 5: Get primer scheme

```
python ./scripts/get_scheme.py \
    --read-file ./results/mpox/fastp/SRR21755837_R1.trim.fastq.gz \
    --scheme-directory ./results/mpox/primerschemes/ \
    yale-mpox/2000/v1.0.0-cladeii
```

### Step 6: Alignment of reads

#### Index the reference genome

```
cp ./results/mpox/primerschemes/yale-mpox/2000/v1.0.0-cladeii/reference.fasta \
    ./results/mpox/bwa/index/
```


```
bwa index \
    -p ./results/mpox/bwa/index/reference \
    ./results/mpox/bwa/index/reference.fasta
```

```
INDEX=$(find -L ./results/mpox/bwa/index -name "*.amb" | sed 's/.amb//')
```


#### Align reads to the reference genome

```
bwa mem \
    -t 4 \
    $INDEX \
    ./results/mpox/fastp/SRR21755837_R1.trim.fastq.gz \
    ./results/mpox/fastp/SRR21755837_R2.trim.fastq.gz | \
    samtools view \
    --threads 2 \
    -bhS \
    -o ./results/mpox/bwa/alignment/SRR21755837.bam - 
```

#### Sort alignment

```
samtools sort \
    --threads 2 \
    -o ./results/mpox/bwa/alignment/SRR21755837.sorted.bam \
    -T ./results/mpox/bwa/alignment/SRR21755837 \
    ./results/mpox/bwa/alignment/SRR21755837.bam
```

#### Index alignment

```
samtools index ./results/mpox/bwa/alignment/SRR21755837.sorted.bam
```

### Step 7: Trim alignments from an amplicon scheme

```
python ./scripts/align_trim.py \
    --normalise 200  \
    ./results/mpox/primerschemes/yale-mpox/2000/v1.0.0-cladeii/primer.bed \
    --paired \
    --no-read-groups \
    --primer-match-threshold 35 \
    --min-mapq 20 \
    --trim-primers \
    --report ./results/mpox/primertrimmed/SRR21755837.alignreport.csv \
    --amp-depth-report ./results/mpox/primertrimmed/SRR21755837.amplicon_depths.tsv \
    < ./results/mpox/bwa/alignment/SRR21755837.sorted.bam 2> ./results/mpox/primertrimmed/SRR21755837.alignreport.er \
    | samtools sort -T SRR21755837 - \
    -o ./results/mpox/primertrimmed/SRR21755837.primertrimmed.rg.sorted.bam \
    && samtools index ./results/mpox/primertrimmed/SRR21755837.primertrimmed.rg.sorted.bam
```

### Step 8: Call variants

```
freebayes \
    -p 1 \
    -f ./results/mpox/bwa/index/reference.fasta \
    -F 0.2 \
    -C 1 \
    --pooled-continuous \
    --min-coverage 10 \
    --gvcf \
    --gvcf-dont-use-chunk true \
    ./results/mpox/primertrimmed/SRR21755837.primertrimmed.rg.sorted.bam \
    > ./results/mpox/freebayes/SRR21755837.gvcf
```

#### Make depth mask, split variants into ambiguous/consensus

```
python ./scripts/process_gvcf.py \
    -d 10 \
    -l 0.25 \
    -u 0.75 \
    -m ./results/mpox/freebayes/SRR21755837.mask.txt \
    -v ./results/mpox/freebayes/SRR21755837.variants.vcf \
    -c ./results/mpox/freebayes/SRR21755837.consensus.vcf \
    ./results/mpox/freebayes/SRR21755837.gvcf
```


#### Normalize variant records into canonical VCF representation

```
for v in "variants" "consensus"; do
    bcftools norm \
        -f ./results/mpox/bwa/index/reference.fasta \
        ./results/mpox/freebayes/SRR21755837.$v.vcf > \
        ./results/mpox/freebayes/SRR21755837.$v.norm.vcf
done
```

### Step 9: Generate consensus genome 
The file is split into a set that should be IUPAC codes and all other bases, using the ConsensusTag in the VCF

```
for vt in "ambiguous" "fixed"; do
    cat ./results/mpox/freebayes/SRR21755837.consensus.norm.vcf | \
    awk -v vartag=ConsensusTag=$vt '$0 ~ /^#/ || $0 ~ vartag' > \
    ./results/mpox/freebayes/SRR21755837.$vt.norm.vcf

    bgzip -f ./results/mpox/freebayes/SRR21755837.$vt.norm.vcf

    tabix -f -p vcf ./results/mpox/freebayes/SRR21755837.$vt.norm.vcf.gz
done
```

#### Apply ambiguous variants first using IUPAC codes. 
This variant set cannot contain indels or the subsequent step will break

```
bcftools consensus \
    -f ./results/mpox/bwa/index/reference.fasta \
    -I ./results/mpox/freebayes/SRR21755837.ambiguous.norm.vcf.gz > \
    ./results/mpox/freebayes/SRR21755837.ambiguous.fa
```

#### Get viral contig name from reference

```
CTG_NAME=$(head -n1 ./results/mpox/bwa/index/reference.fasta | sed 's/>//')
```


#### Apply remaning variants, including indels

```
bcftools consensus \
    -f ./results/mpox/freebayes/SRR21755837.ambiguous.fa \
    -m ./results/mpox/freebayes/SRR21755837.mask.txt \
    ./results/mpox/freebayes/SRR21755837.fixed.norm.vcf.gz | \
    sed s/$CTG_NAME/SRR21755837/ > \
    ./results/mpox/freebayes/SRR21755837.consensus.fa
```

#### Coverage metrics and visualization with IGV

### Step 10: Squirrel - Some QUIck Reconstruction to Resolve Evolutionary Links
The MPXV genome is pretty challenging to work with and do reliable phylogenetics
on. It is large (~200kb), has tracts of low complexity and repetitive regions,
and has large deletions, which can lead to difficulties producing a reliable
alignment. With squirrel, we provide a rapid way of producing reliable
alignments for MPXV and also enable maximum-likelihood phylogenetics pipeline
tree estimation.

```
export XDG_CACHE_HOME=$PWD/.cache
```

#### Add additional sequences retrieved from Pathoplexus

```
wget -c https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/mpoxv/mpox_nuc_DRC_2024.fasta \
-P ./results/mpox/pathoplexus-data/
```

```
wget -c https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/mpoxv/mpox_metadata_DRC_2024.txt \
-P ./results/mpox/pathoplexus-data/
```

# Concatenate consensus fasta files

```
cat \
    ./results/mpox/pathoplexus-data/*.fasta \
    ./results/mpox/freebayes/SRR21755837.consensus.fa \
    > ./results/mpox/all_consensus/mpxv_all_consensus.fasta
```

Enrichment of APOBEC3-mutations in the MPXV population are a signature of
sustained human-to-human transmission. Identifying APOBEC3-like mutations in
MPXV genomes from samples in a new outbreak can be a piece of evidence to
support sustained human transmission of mpox. Squirrel can run an
APOBEC3-reconstruction and map these mutations onto the phylogeny.

Squirrel maps each query genome in the input file against a reference genome specific to each clade using minimap2. Using gofasta, the mapping file is then converted into a multiple sequence alignment.

For Clade II, the reference used is NC_063383 and for Clade I, we use NC_003310.
This means that all coordinates within an alignment will be relative to these
references. A benefit of this is that within a clade, alignment files and be
combined without having to recalculate the alignment.

```
squirrel \
    results/mpox/all_consensus/mpxv_all_consensus.fasta \
    --no-mask \
    --seq-qc \
    -o squirrel \
    --outfile all_consensus.aln.fasta \
    --tempdir results/mpox/tmp/squirrel/ \
    -t 2 \
    --run-phylo \
    --run-apobec3-phylo \
    --outgroups KJ642617,KJ642615 \
    --clade split
```
