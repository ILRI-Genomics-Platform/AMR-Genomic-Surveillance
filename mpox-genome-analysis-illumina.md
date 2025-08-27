
# *Monkeypox virus* analysis workflow using Illumina data
---  

###### **_Trainers_**: [John Juma](https://github.com/ajodeh-juma), [Kennedy Mwangi](https://github.com/wanjauk), [Ouso Daniel](https://github.com/ousodaniel) & [Gilbert Kibet](https://github.com/kibet-gilbert)

---

## Table of Contents
  <!-- - [Overview](#overview)
  - [Learning Objectives](#learning-objectives) -->
  - [Background](#background)
  - [Epidemiology in Africa](#epidemiology-in-africa)
  - [About the sample](#about-the-sample)
  - [Set-up](#set-up)
      - [Workshop Environment](#workshop-environment)
      - [Logging into the HPC](#logging-into-the-hpc)
      - [Compute Node](#compute-node)
      - [Project organisation](#project-organisation)
  - [Bioinformatics Analysis](#bioinformatics-analysis)
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
      - [Step 8: Call variants](#step-8-call-variants)
        - [Make depth mask, split variants into ambiguous and consensus](#make-depth-mask-split-variants-into-ambiguous-and-consensus)
        - [Normalize variant records into canonical VCF representation](#normalize-variant-records-into-canonical-vcf-representation)
      - [Step 9: Generate consensus genome](#step-9-generate-consensus-genome)
        - [Apply ambiguous variants first using IUPAC codes](#apply-ambiguous-variants-first-using-iupac-codes)
        - [Get viral contig name from reference](#get-viral-contig-name-from-reference)
        - [Apply remaning variants, including
          indels](#apply-remaning-variants-including-indels)
        - [Coverage metrics and visualization with IGV](#coverage-metrics-and-visualization-with-igv)
      - [Step 10: Squirrel Some QUIck Reconstruction to Resolve Evolutionary Links](#step-10-squirrel-some-quick-reconstruction-to-resolve-evolutionary-links)
        - [Add additional sequences retrieved from
          Pathoplexus](#add-additional-sequences-retrieved-from-pathoplexus)
        - [Fetch genome sequences](#fetch-genome-sequences)
      

# Background

Monkeypox virus (MPXV) was first isolated in Denmark in the late 1950s from a
colony of laboratory monkeys from Singapore that were going to be used for polio
virus research. During the following decade, additional outbreaks of mpox were
seen in laboratory animals in the United States as well as zoo animals in
Rotterdam. MPXV was first identified as a cause of disease in humans in the
1970s in what is now the Democratic Republic of the Congo (DRC).


In November 2022, the World Health Organization, who is responsible for naming and renaming of diseases under the International Classification of Diseases (ICD), changed the name of the disease referred to as “monkeypox” to “mpox” [1]. This change was made to follow current best practices of not naming diseases after animals or geographic locations and to reduce any stigma that could be associated with the original name.

Mpox is an infectious disease caused by the monkeypox virus (MPXV). There are
two known clades of MPXV: clade I (previously called the Congo Basin clade),
which includes subclades Ia and Ib; and clade II (previously called the West
Africa clade), which includes subclades IIa and clade IIb. Historically mpox was
primarily characterized by spread from animals to humans, but in recent years,
more and more human-to-human transmission has occurred. It can spread among
humans through: 
- direct close physical contact with an infected person, including sexual contact
- indirect contact, that is, through contact with contaminated materials
- respiratory contact, through infectious respiratory particles
- mother-to-child transmission (vertical transmission).

The existence of subclades of clade I was detected during the outbreak in the
Democratic Republic of the Congo (DRC) that started in 2023. 

# Epidemiology in Africa

Outbreak in Central and East Africa starting in 2023 — On August 14, 2024, the
World Health Organization (WHO) declared for a second time that mpox was a
global public health emergency[REF]. From January 1, 2024 to August 18, 2024, a
total of 21,786 laboratory-confirmed mpox cases leading to 607 deaths were
reported from 12 African countries, with approximately 90 percent of the cases
occurring in the DRC [REF].


# About the Sample
**Project accession:** [PRJNA885473](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA885473)


# Set-up  
## Workshop Environment
The workshop is mainly based on working in a HPC environment. However, local
set-ups will be provisioned, although some steps of the workflows may be
impossible to execute on laptops. For access to the HPC, you will use your
personal computers to log into the ILRI HPC cluster, which operates on a Linux
operating system. If using a non-Unix operating system you will require a
program that enables interfacing with the Linux HPC environment.


## Logging into the HPC  
To log in to the HPC, you use the provided details (username and password) with the below command. The username follows the parttern `Bio4InfoXX`, where `XX` is a number.
```
ssh <user_name>@hpc.ilri.cgiar.org
```
- Replace `<user_name>` in the command with the provided username and execute (press enter). 
- Next, you will be prompted for the password.  
 >***Note:*** The password will not be visible to you as you type--just a little faith.

## Compute Node
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

## Project Organisation  
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

2. *Change directory into the created directory*

    ```
    cd /var/scratch/$USER/ACDC_AMR2025
    ```

>**Note** Once inside the `hpc`, all instances of ```$USER``` will be equivalent to the `hpc` username that you were assigned. Your username, by default, is stored in a variable called `USER`. By using it, you will not have to type your username, rather, your shell will automatically retrieve your username which is the value stored in the `USER` variable. The `$` (dollar) prefix is used to retrieve the value of that variable.

3. *Create project sub-directories:*

    ```
    mkdir -p \
    results/mpox/{fastqc,fastp,hostile,trim_galore,primerschemes,bwa/{index,alignment},tmp/squirrel,primertrimmed,primer-trimmed,freebayes,data/{ncbi,pathoplexus,all-consensus},squirrel}
    ```

3. *Create symblic links to the required resources*

    ```
    ln -sf /var/scratch/global/jjuma/ACDC_AMR2025/[dpsr]* . 
    ```


4. *Retrieve files from the GitHub repository*

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

# Bioinformatics Analysis

# Step 1: Load required modules

```
module load hostile/2.0.0
module load fastp/0.22.0
module load seqtk/1.3
module load bwa/0.7.17
module load freebayes/1.3.4
module unload bcftools/1.17
module load bcftools/1.8
module load squirrel/1.1.2
```

# Step 2: Create a virtual environment and install dependencies

```
python3 -m venv ./py3env
```

```
source ./py3env/bin/activate
```

```
python3 -m pip install -r ./utils/artic-requirements.txt
```

# Step 3: Remove human reads

```
hostile clean \
    --fastq1 ./data/mpox/illumina/SRR21755837_1.fastq.gz \
    --fastq2 ./data/mpox/illumina/SRR21755837_2.fastq.gz \
    --output ./results/mpox/hostile/ \
    --threads 4 \
    --index ./databases/hostile/human-t2t-hla
```

# Step 4: Trim adapter sequences

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
 
# Step 5: Get primer scheme

```
cp ./primer_scheme/reference.fasta \
    ./primer_scheme/primer.bed \
    ./results/mpox/primerschemes/
```

<!-- ```
python ./scripts/get_scheme.py \
    --read-file ./results/mpox/fastp/SRR21755837_R1.trim.fastq.gz \
    --scheme-directory ./results/mpox/primerschemes/ \
    yale-mpox/2000/v1.0.0-cladeii
``` -->

# Step 6: Alignment of reads

## Index the reference genome

<!-- ```
cp ./results/mpox/primerschemes/yale-mpox/2000/v1.0.0-cladeii/reference.fasta \
    ./results/mpox/bwa/index/
``` -->


```
bwa index \
    -p ./results/mpox/bwa/index/reference \
    ./results/mpox/primerschemes/reference.fasta
```

```
INDEX=$(find -L ./results/mpox/bwa/index -name "*.amb" | sed 's/.amb//')
```


## Align reads to the reference genome

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

## Sort alignment

```
samtools sort \
    --threads 2 \
    -o ./results/mpox/bwa/alignment/SRR21755837.sorted.bam \
    -T ./results/mpox/bwa/alignment/SRR21755837 \
    ./results/mpox/bwa/alignment/SRR21755837.bam
```

## Index alignment

```
samtools index ./results/mpox/bwa/alignment/SRR21755837.sorted.bam
```

# Step 7: Trim alignments from an amplicon scheme

1. Given a reference position and a direction of travel, walk out and find the nearest primer site.
2. Soft mask an alignment to fit within primer start/end sites
    - softmask the alignment if left primer start/end inside alignment
    - softmask the alignment if right primer start/end inside alignment
3. Handle the alignment segment including
    - filter out unmapped and supplementary alignment segments
    - locate the nearest primers to this alignment segment
    - check if primers are correctly paired and then assign read group
    - get the amplicon number
    - report with this alignment segment + primer details

4. Generate a dictionary of amplicons from a primer scheme list (generated by
   vcftagprimersites/read_bed_file)
5. Normalise the depth of the trimmed segments to a given value. Perform
   per-amplicon normalisation using numpy vector maths to determine whether the
   segment in question would take the depth closer to the desired depth across
   the amplicon.
6. Generate read pairs in a BAM file or within a region string. Reads are added
   to read_dict until a pair is found.
7. Filter and soft mask an alignment file so that the alignment boundaries match
   the primer start and end sites.Based on the most likely primer position,
   based on the alignment coordinates.


```
python ./scripts/align_trim.py \
    --normalise 200  \
    ./results/mpox/primerschemes/primer.bed \
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

# Step 8: Call variants

freebayes is a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment.

freebayes is haplotype-based, in the sense that it calls variants based on the
literal sequences of reads aligned to a particular target, not their precise
alignment. This model is a straightforward generalization of previous ones (e.g.
PolyBayes, samtools, GATK) which detect or report variants based on alignments.
This method avoids one of the core problems with alignment-based variant
detection--- that identical sequences may have multiple possible alignments.


```
freebayes \
    -p 1 \
    -f ./results/mpox/primerschemes/reference.fasta \
    -F 0.2 \
    -C 1 \
    --pooled-continuous \
    --min-coverage 10 \
    --gvcf \
    --gvcf-dont-use-chunk true \
    ./results/mpox/primertrimmed/SRR21755837.primertrimmed.rg.sorted.bam \
    > ./results/mpox/freebayes/SRR21755837.gvcf
```

`--pooled-continuous`: Output all alleles which pass input filters/thresholds,
regardless of genotyping outcome or model.

`--min-coverage`: -Require at least this coverage to process a site.

## Make depth mask, split variants into ambiguous and consensus

Process a .gvcf file to create a file of consensus variants, low-frequency
variants and a coverage mask

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

`-d`: Mask reference positions with depth less than this threshold

`-l`: Variants with frequency less than -l will be discarded

`-u`: Substitution variants with frequency less than -u will be encoded with
IUPAC ambiguity codes

`-v`: The output file name for variants (non-reference gVCF records)

`-c`: The output file name for variants that will be applied to generate the
consensus sequence



## Normalize variant records into canonical VCF representation
Realign and normalize indels; check if REF alleles match the reference;
split multiallelic sites into multiple rows; recover multiallelics from
multiple rows.

```
for v in "variants" "consensus"; do
    echo -e "normalising variants in: $v"
    bcftools norm \
        -f ./results/mpox/primerschemes/reference.fasta \
        ./results/mpox/freebayes/SRR21755837.$v.vcf > \
        ./results/mpox/freebayes/SRR21755837.$v.norm.vcf
done
```

total/split/joined/realigned/skipped

# Step 9: Generate consensus genome 
The consensus VCF file is split into a set that should be IUPAC codes and all
other bases, using the ConsensusTag in the VCF

```
for vt in "ambiguous" "fixed"; do
    echo "Splitting on ConsensusTag on: $vt"
    cat ./results/mpox/freebayes/SRR21755837.consensus.norm.vcf | \
    awk -v vartag=ConsensusTag=$vt '$0 ~ /^#/ || $0 ~ vartag' > \
    ./results/mpox/freebayes/SRR21755837.$vt.norm.vcf

    bgzip -f ./results/mpox/freebayes/SRR21755837.$vt.norm.vcf

    tabix -f -p vcf ./results/mpox/freebayes/SRR21755837.$vt.norm.vcf.gz
done
```

## Apply ambiguous variants first using IUPAC codes. 
This variant set cannot contain indels or the subsequent step will break

```
bcftools consensus \
    -f ./results/mpox/primerschemes/reference.fasta \
    -I ./results/mpox/freebayes/SRR21755837.ambiguous.norm.vcf.gz > \
    ./results/mpox/freebayes/SRR21755837.ambiguous.fa
```

## Get viral contig name from reference

```
CTG_NAME=$(head -n1 ./results/mpox/primerschemes/reference.fasta | sed 's/>//')
```


## Apply remaning variants, including indels

```
bcftools consensus \
    -f ./results/mpox/freebayes/SRR21755837.ambiguous.fa \
    -m ./results/mpox/freebayes/SRR21755837.mask.txt \
    ./results/mpox/freebayes/SRR21755837.fixed.norm.vcf.gz | \
    sed s/$CTG_NAME/SRR21755837/ > \
    ./results/mpox/freebayes/SRR21755837.consensus.fa
```

## Coverage metrics and visualization with IGV

Transfer the `reference.fa` and `reference.fa.fai` and the sorted primertimmed
alignment files to hpc

```
cp ./results/mpox/primerschemes/reference.fa* ~/
```

```
cp ./results/mpox/primertrimmed/SRR21755837.primertrimmed.* ~/
```

Now transfer the files to your local computer

```
rsync -avP --partial <user_name>@hpc.ilri.cgiar.org:~/reference.fa* ~/
```

```
rsync -avP --partial <user_name>@hpc.ilri.cgiar.org:~/SRR21755837.primertrimmed.* ~/
```

Open the IGV web app and load the reference genome and the alignment

Zoom into the region `MT903345:7,766-7,804`. 
What SNP do you observe? 
What is the position of the SNP?
What is the alternate allele, what is the reference allele?
How many alignments span this region?

Now enter into the location `MT903345:12,835-12,874`
What SNP do you observe? 
What is the position of the SNP?
What is the alternate allele, what is the reference allele?
How many alignments span this region?


# Step 10: Squirrel Some QUIck Reconstruction to Resolve Evolutionary Links
The MPXV genome is pretty challenging to work with and do reliable phylogenetics
on. It is large (~200kb), has tracts of low complexity and repetitive regions,
and has large deletions, which can lead to difficulties producing a reliable
alignment. With squirrel, we provide a rapid way of producing reliable
alignments for MPXV and also enable maximum-likelihood phylogenetics pipeline
tree estimation.

<!-- ```
export XDG_CACHE_HOME=$PWD/.cache
``` -->

## Add additional sequences retrieved from Pathoplexus

```
wget -c https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/mpoxv/mpox_nuc_DRC_2024.fasta \
-P ./results/mpox/data/pathoplexus
```


```
wget -c https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/mpoxv/mpox_metadata_DRC_2024.txt \
-P ./results/mpox/data/pathoplexus/
```


```
wget -c https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/mpoxv/MPOXV-accessions.txt \
-P ./results/mpox/data/ncbi/
```


## Fetch genome sequences

**Project paper:** [Publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC11484285/pdf/eurosurv-29-38-1.pdf)
We will use 12 MPOXV genomes sequenced from the DRC. We can either retrieve these
sequences from the NCBI or Pathoplexus.


```
python ./scripts/fetch-genomes.py \
    ./results/mpox/data/ncbi/MPOXV-accessions.txt \
    ./results/mpox/data/ncbi/
```

## Concatenate consensus fasta files

```
cat \
    ./results/mpox/data/ncbi/KJ6* \
    ./results/mpox/data/pathoplexus/*.fasta \
    ./results/mpox/freebayes/SRR21755837.consensus.fa \
    > ./results/mpox/data/all-consensus/all_consensus.fasta
```

```
sed 's/\.1.*//g' ./results/mpox/data/all-consensus/all_consensus.fasta \
    > ./results/mpox/data/all-consensus/mpxv_all_consensus.fasta
```

## Run Squirrel

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


<!-- #### Create a conda environemnt for squirrel

```
module purge
```

```
deactivate
```

```
mamba create -c bioconda -c conda-forge -n squirrel -y squirrel
```

#### Activate the conda environment
```
conda activate squirrel
``` -->


```
squirrel \
    results/mpox/data/all-consensus/mpxv_all_consensus.fasta \
    --no-mask \
    --seq-qc \
    --outdir results/mpox/squirrel \
    --outfile all_consensus.aln.fasta \
    --threads 1 \
    --run-phylo \
    --run-apobec3-phylo \
    --interactive-tree \
    --outgroups KJ642617,KJ642615,KJ642616 \
    --clade cladei
```

## Output files

Column | Description
--- | ---
`.aln.fasta` | The alignment file, with alignment scaffolded against a clade-specific reference. By default one of the ITR regions and a curated set of problematic regions is masked as Ns.
`.aln.tree` | The output maximum likelihood tree file from IQTREE2 with Node labels that correspond to the reconstruction Node labels. This tree can be viewed in various tree viewers, for example FigTree.
`.aln.tree.state` / `.aln.tree.state_differences.csv` | The output ancestral state reconstruction file from IQTREE2 and the compiled list of unambiguously variable sites from squirrel.
`.aln.tree.branch_snps.reconstruction.csv` | A report of individual site changes mapped to specific branches and their dinucleotide context.
`.aln.tree.amino_acid.reconstruction.csv` | A report of each mutation that occurs across the phylogeny, their location, dinucleotide context, APOBEC3 status, which gene they're present in, codon position, amino acid change and a prediction of how extreme that amino acid change is with Grantham score.
`.aln.tree.png` / `sequences.aln.tree.svg` |  Visualisation of reconstructed tree showing whether mutations are consistent with APOBEC3 editing or not.
`.aln.report.html` | Summary report of analysis run.


## Copy output files to `home` directory

```
rsync -avP --partial ./results/mpox/squirrel ~/ --exclude="*.fasta" --exclude="*.state"
```

## Copy output files from `home` directory to local machine/compute

```
rsync -avP --partial <user_name>@hpc.ilri.cgiar.org:~/squirrel ~/ --exclude="*.fasta" --exclude="*.state"
```
<!-- 
# Temporal signal

```
tail -n+2 ./results/mpox/data/pathoplexus/mpox_metadata_DRC_2024.txt | cut -f1,14 > ./results/mpox/data/pathoplexus/dates.txt
``` -->