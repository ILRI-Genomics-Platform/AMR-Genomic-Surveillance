
```
module load resfinder/4.6.0
```

## Get help
```
python -m resfinder -h
```

## Run - ecoli
```
python -m resfinder \\
    -ifa data/fasta/ecoli/MSAK01.fasta \\
    -o data/resfinder/ecoli \\
    -s ecoli \\
    --min_cov 0.6 \\
    --threshold 0.9 \\
    --min_cov_point 0.6 \\
    --threshold_point 0.9 \\
    --ignore_stop_codons \\
    --ignore_indels \\
    --acquired --point
```

## Run - klebs
```
python -m resfinder \\
    -ifa data/fasta/klebs/MSAK01.fasta \\
    -o data/resfinder/klebs \\
    -s klebsiella \\
    --min_cov 0.6 \\
    --threshold 0.9 \\
    --min_cov_point 0.6 \\
    --threshold_point 0.9 \\
    --ignore_stop_codons \\
    --ignore_indels \\
    --acquired --point
```
