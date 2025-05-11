
# *Monkeypox virus* analysis workflow using Illumina data
---  

###### **_Trainers_**: [John Juma](https://github.com/ajodeh-juma), [Kennedy Mwangi](https://github.com/wanjauk), [Ouso Daniel](https://github.com/ousodaniel) & [Gilbert Kibet](https://github.com/kibet-gilbert)

---

# Set up directories
Before starting the analysis, ensure that you are logged into the HPC, create an interactive session on the assigned compute node, and change directory to the project folder which is `ACDC_AMR2025`.

```
mkdir -p \
results/mpox/{fastqc,fastp,hostile,trim_galore,primerschemes,bwa/{index,alignment},primertrimmed,primer-trimmed,freebayes}

ln -sf /var/scratch/global/jjuma/ACDC_AMR2025/[dpsr]* . 

cp /var/scratch/global/jjuma/ACDC_AMR2025/artic-requirements.txt . 
```

# Load modules
```
module load hostile/2.0.0
module load fastp/0.22.0
module load seqtk/1.3
module load bwa/0.7.17
module load freebayes/1.3.4
module unload bcftools/1.17
module unload bcftools/1.13
```

## Create a virtual and install dependencies

```
python3 -m venv $BASEDIR/py3env
```

```
source $BASEDIR/py3env/bin/activate
```

```
python3 -m pip install -r ./artic-requirements.txt
```

# Remove human reads

```
hostile clean \
    --fastq1 ./data/mpox/illumina/SRR21755837_1.fastq.gz \
    --fastq2 ./data/mpox/illumina/SRR21755837_2.fastq.gz \
    --output ./results/mpox/hostile/ \
    --threads 4 \
    --index ./databases/hostile/human-t2t-hla
```

<<<<<<< HEAD
# trim adapters using trim galore
=======
# Trim adapters using trim galore
>Note cutadapt which is a dependency of trim galore is not properly configured on hpc,
we opt to use fastp instead
>>>>>>> 8a140af98364322a064007f5effa221cc82b859e

```
trim_galore \
    --cores 2 \
    --output_dir ./results/mpox/trim_galore \
    --paired \
    ./results/mpox/hostile/SRR21755837_1.clean_1.fastq.gz \
    ./results/mpox/hostile/SRR21755837_2.clean_2.fastq.gz
```
<<<<<<< HEAD
=======

# Trim adapters using fastp

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
>>>>>>> 8a140af98364322a064007f5effa221cc82b859e
 
## Get scheme

```
python ./scripts/get_scheme.py \
    --read-file ./results/mpox/fastp/SRR21755837_R1.trim.fastq.gz \
    --scheme-directory ./results/mpox/primerschemes/ \
    yale-mpox/2000/v1.0.0-cladeii
```

# Alignment

## Index the reference genome

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


## Align

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

## Sort

```
samtools sort \
    --threads 2 \
    -o ./results/mpox/bwa/alignment/SRR21755837.sorted.bam \
    -T ./results/mpox/bwa/alignment/SRR21755837 \
    ./results/mpox/bwa/alignment/SRR21755837.bam
```


## Index

```
samtools index ./results/mpox/bwa/alignment/SRR21755837.sorted.bam
```

## Trim alignments from an amplicon scheme

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

## Call variants

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

## Make depth mask, split variants into ambiguous/consensus

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


## Normalize variant records into canonical VCF representation

```
for v in "variants" "consensus"; do
    bcftools norm \
        -f ./results/mpox/bwa/index/reference.fasta \
        ./results/mpox/freebayes/SRR21755837.$v.vcf > \
        ./results/mpox/freebayes/SRR21755837.$v.norm.vcf
done
```

## Split the consensus sites file 
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

## Apply ambiguous variants first using IUPAC codes. 
This variant set cannot contain indels or the subsequent step will break

```
bcftools consensus \
    -f ./results/mpox/bwa/index/reference.fasta \
    -I ./results/mpox/freebayes/SRR21755837.ambiguous.norm.vcf.gz > \
    ./results/mpox/freebayes/SRR21755837.ambiguous.fa
```

## Get viral contig name from reference

```
CTG_NAME=$(head -n1 ./results/mpox/bwa/index/reference.fasta | sed 's/>//')
```


## Apply remaninng variants, including indels

```
bcftools consensus \
    -f ./results/mpox/freebayes/SRR21755837.ambiguous.fa \
    -m ./results/mpox/freebayes/SRR21755837.mask.txt \
    ./results/mpox/freebayes/SRR21755837.fixed.norm.vcf.gz | \
    sed s/$CTG_NAME/SRR21755837/ > \
    ./results/mpox/freebayes/SRR21755837.consensus.fa
```

## Coverage metrics and visualization with IGV

## Squirrel

```
export XDG_CACHE_HOME=$PWD/.cache
```
