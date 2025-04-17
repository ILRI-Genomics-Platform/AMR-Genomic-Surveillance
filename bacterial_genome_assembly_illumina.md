
# QC
fastp \
    --in1 SRR25008769_R1.fastq.gz \
    --in2 SRR25008769_R2.fastq.gz \
    --out1 filt-r1.fq \
    --out2 filt-r2.fq \
    --detect_adapter_for_pe \
    --thread 2 \
    --json results/summary/SRR25008769.fastp.json \
    --html results/summary/SRR25008769.fastp.html \ --cut_mean_quality 20 \
    --qualified_quality_phred 25 \
    --unqualified_percent_limit 40 \
    --length_required 20 2> SRR25008769-fastp.log



# assembly

shovill \
    --R1 SRR25008770_R1.fastq.gz \
    --R2 SRR25008770_R2.fastq.gz \
    --gsize 5249449 \
    --outdir results \
    --assembler skesa \
    --minlen 500 \
    --mincov 2 \
    --force \
    --keepfiles \
    --depth 0 \
    --noreadcorr \
    --namefmt "SRR25008770_%05d" \
    --cpus 2 \
    --ram 7


# Annotation
export PROKKA_DBDIR=$(echo "$(which prokka | sed "s=/prokka==")/../db")
env
mkdir tmp_prokka/
TMPDIR=tmp_prokka/ bactopia-prokka \
    --evalue 1e-09 --coverage 80 --centre Bactopia \
    --cpus 2 \
    --prefix SRR25008770 \
    --locustag SRR25008770 \
    --proteins proteins.faa \
    SRR25008770.fna

rm -rf tmp_prokka/


AMRFINDER_DB=$(find amrfinderplus.tar.gz/ -name "AMR.LIB" | sed 's=AMR.LIB==')

# Full AMRFinderPlus search combining results
amrfinder \
    --nucleotide SRR25008770.fna \
    --protein SRR25008770.faa \
    --gff SRR25008770.gff \
    --annotation_format prokka \
    --plus \
    --ident_min -1 \
    --coverage_min 0.5 \
    --translation_table 11 \
    --database $AMRFINDER_DB \
    --threads 2 \
    --name SRR25008770 > SRR25008770.tsv

# Multilocus sequence typing
MLST_DB=$(find mlst.tar.gz/ -name "mlst.fa" | sed 's=blast/mlst.fa==')
mlst \
    --threads 2 \
    --blastdb $MLST_DB/blast/mlst.fa \
    --datadir $MLST_DB/pubmlst \
    --scheme ecoli_achtman_4 \
    --minid 95 \
    --mincov 10 \
    --minscore 50 \
    SRR25008770.fna.gz \
    > SRR25008770.tsv


