
# Introduction to High-Performance Computing with SLURM
---

## HPC Concepts

**High-Performance Computing (HPC)**

High-Performance Computing refers to using computing clusters for computationally intensive tasks. Key benefits for bioinformatics include:

- Processing large datasets that exceed desktop capabilities
- Running analyses in parallel: several CPUs per Job and several Jobs
- Access to more memory, storage, and computing power
- Specialized hardware for specific tasks

**Cluster Architecture**

A typical HPC cluster consists of:

- **Login nodes**: Where users connect and manage jobs
- **Compute nodes**: Where jobs actually run
- **Storage systems**: For data and results
- **Interconnect**: High-speed network connecting nodes
- **Job scheduler**: Software that manages job execution

**Resource Management Concepts**

- **Jobs**: Computational tasks submitted to the cluster
- **Resources**: CPUs, memory, time, GPUs, etc.
- **Queues/Partitions**: Job categories with different priorities and resources
- **Allocation**: Amount of resources available to a user/group
- **Walltime**: Maximum execution time for a job

## SLURM Basics

**What is SLURM?**

SLURM (Simple Linux Utility for Resource Management) is a workload manager for Linux clusters. It provides three key functions:

1. Allocate resources to users
2. Provide a framework for starting/monitoring jobs
3. Manage a queue of pending jobs

**Key SLURM Commands**

- `sinfo`: View information about compute nodes and partitions
  ```bash
  # Show cluster status
  sinfo
  
  # Show detailed node information
  sinfo -N -l
  ```

- `squeue`: View information about jobs in the queue
  ```bash
  # Show all jobs
  squeue
  
  # Show only your jobs
  squeue -u $USER
  ```

- `sbatch`: Submit a batch job script
  ```bash
  sbatch my_job.sh
  ```

- `scancel`: Cancel a job
  ```bash
  # Cancel a specific job
  scancel 12345
  
  # Cancel all your jobs
  scancel -u $USER
  ```

## Creating SLURM Job Scripts

 - Before we proceed, first [log in to the HPC](https://github.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/tree/main#logging-into-the-hpc) as explained earlier.  

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

We will be conducting our a anlysis from a directory in the `scratch` space of the HPC. It is a temporary storage area designed for high-speed data access and short-term file storage. It is cleaned continually (say every night) of files older than X time (90 days).
```
mkdir -p /var/scratch/$USER
cd /var/scratch/$USER
```

 - Then let us copy some test data.  
 
Our test data is in the path `/var/scratch/global/gkibet/ACDC_AMR2025/data/test_data`. We will copy the folder to current directory.

```bash
pwd # Ensure you are in /var/scratch/$USER/

# Copy data
cp -rf /var/scratch/global/gkibet/ACDC_AMR2025/data/test_data/ ./

# Check that data is copied
ls test_data

# List the contents of the directory to make sure all files are there:  
ls test_data/*/
```

**Basic SLURM Script Structure**

To submit a job to the `HPC` via `slurm` job scheduler you may use a script. The script will list the parameters that slurm needd to process that job and the commands it will execute. This script is reffered to as an **sbatch script**. It may be given a suffix like `.sbatch` but is ideally a bash executable script.

We can break down an sbatch script into the following segments: 

 - shebang  

```bash
#!/bin/bash
```
 - sbatch parameters - set by starting a line with `#SBATCH` followed by the option and argument: `#SBATCH --OPTION=ARGUMENT`.

```bash
#SBATCH --job-name=myanalysis      # Job name
#SBATCH --output=myanalysis_%j.out # Output file (%j expands to job ID)
#SBATCH --error=myanalysis_%j.err  # Error file
#SBATCH --time=01:00:00            # Wall time (HH:MM:SS)
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --ntasks=1                 # Number of tasks
#SBATCH --cpus-per-task=4          # Number of CPUs per task
#SBATCH --mem=8G                   # Memory per node
#SBATCH --partition=normal         # Partition/queue name
```

 - Loading specific software modules and versions

```bash
# Load any required modules
module load samtools/1.17 bwa/0.7.19
```

 - Preliminary bash commands - helps report on errors or basic settings as set in sbatch segment above.

```bash
# Run commands
echo "Job started at $(date)"
echo "Running on host: $(hostname)"
echo "Working directory: $(pwd)"
```

 - Actual bioinformatics analysis commands

```bash
# Use bwa to index a reference genome
bwa index test_data/reference/reference.fa

# Example bioinformatics commands
bwa mem -t $SLURM_CPUS_PER_TASK \
	./test_data/reference/reference.fa \
	./test_data/fastq/read1.fq \
	./test_data/fastq/read2.fq > output.sam

samtools sort -@ $SLURM_CPUS_PER_TASK -o output.bam output.sam
samtools index output.bam
```
 - Post analysis bash commands e.g reporting results or copying it to more permanent storage
```bash
echo "Job finished at $(date)"
```

<details>
    <summary>
        Click to toggle contents of 
        <b style='color: blue'> Actual script
        </b>
    </summary>
    When put together the sbatch script above looks like:

```bash
#!/bin/bash
#SBATCH --job-name=myanalysis      # Job name
#SBATCH --output=myanalysis_%j.out # Output file (%j expands to job ID)
#SBATCH --error=myanalysis_%j.err  # Error file
#SBATCH --time=01:00:00            # Wall time (HH:MM:SS)
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --ntasks=1                 # Number of tasks
#SBATCH --cpus-per-task=4          # Number of CPUs per task
#SBATCH --mem=8G                   # Memory per node
#SBATCH --partition=normal         # Partition/queue name

# Load any required modules
module load samtools/1.17 bwa/0.7.19

# Run commands
echo "Job started at $(date)"
echo "Running on host: $(hostname)"
echo "Working directory: $(pwd)"

# Use bwa to index a reference genome
bwa index test_data/reference/reference.fa

# Example bioinformatics commands
bwa mem -t $SLURM_CPUS_PER_TASK \
	./test_data/reference/reference.fa \
	./test_data/fastq/read1.fq \
	./test_data/fastq/read2.fq > output.sam

samtools sort -@ $SLURM_CPUS_PER_TASK -o output.bam output.sam
samtools index output.bam

echo "Job finished at $(date)"
```

</details>


<details>
    <summary>
        Click to toggle contents of 
        <b style='color: blue'> SLURM Directives
        </b>
    </summary>
    Common SLURM directives for bioinformatics jobs:

| Directive | Description | Example |
|-----------|-------------|---------|
| `--job-name` | Name of the job | `--job-name=alignment` |
| `--output` | Standard output file | `--output=align_%j.out` |
| `--error` | Standard error file | `--error=align_%j.err` |
| `--time` | Maximum wall time | `--time=12:00:00` |
| `--nodes` | Number of nodes | `--nodes=1` |
| `--ntasks` | Number of tasks | `--ntasks=1` |
| `--cpus-per-task` | CPUs per task | `--cpus-per-task=8` |
| `--mem` | Memory per node | `--mem=32G` |
| `--mem-per-cpu` | Memory per CPU | `--mem-per-cpu=4G` |
| `--partition` | Partition to use | `--partition=compute` |
| `--account` | Account to charge | `--account=bioinformatics` |
| `--mail-user` | Email for notifications | `--mail-user=user@example.com` |
| `--mail-type` | When to send email | `--mail-type=END,FAIL` |
| `--array` | Job array indices | `--array=1-10` |
| `--dependency` | Job dependencies | `--dependency=afterok:12345` |

</details>

**Environment Variables**

SLURM sets environment variables that can be used in your script:

```bash
echo "Job ID: $SLURM_JOB_ID"
echo "Number of nodes: $SLURM_JOB_NUM_NODES"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Node list: $SLURM_JOB_NODELIST"
```

---
### Bioinformatics Workflows with SLURM
---

#### Running Bioinformatics Tools with SLURM

Now let us practice some scripting. In this case we will create a script that will submit a "**job**" through `slurm` to the HPC.

**Problem statement:** We have some reads and we need to align them to a reference genome and generate some summary statistics. To break it down further, our script will have the following features:  

 - Start with a shebang.  
 - Set some `slurm` `sbatch` parameters.  
 - Load the needed modules.
 - Set the paths to the reference genome and the reads as variables.
 - Make sure you have created an output directory.
 - Record the start time for the analysis.
 - Peform the actual analysis: **index the genome**, **align the reads**, **index the BAM file**, **calculate the analysis statistics**.   
 - Record the end time of the analysis and calculate how long it took and report.

**Reference genome:** `test_data/reference/genome.fa`  
**Input Reads:** `test_data/fastq/sample1_R1.fastq.gz` and `test_data/fastq/sample1_R2.fastq.gz`  
**Tools:** `samtools/1.17` and `bwa/0.7.19`  

**Example: DNA Sequence Alignment**

```bash
#!/bin/bash
#SBATCH --job-name=bwa_align
#SBATCH --output=bwa_align_%j.log
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# Load modules
module load samtools/1.17 bwa/0.7.19

# Set variables
REFERENCE="./test_data/reference/genome.fa"
SAMPLE_ID="sample1"
READ1="./test_data/fastq/${SAMPLE_ID}_R1.fastq.gz"
READ2="./test_data/fastq/${SAMPLE_ID}_R2.fastq.gz"
OUTPUT_DIR="./test_data/output"

# Create output directory
mkdir -p $OUTPUT_DIR

# Start timing
start_time=$(date +%s)

# Indexing the reference 
echo "Indexing reference genome in readiness for alignment"
bwa index $REFERENCE

# Run alignment
echo "Starting alignment for $SAMPLE_ID"
bwa mem -t $SLURM_CPUS_PER_TASK $REFERENCE $READ1 $READ2 | \
    samtools sort -@ $SLURM_CPUS_PER_TASK -o $OUTPUT_DIR/${SAMPLE_ID}.bam -

# Index BAM file
samtools index $OUTPUT_DIR/${SAMPLE_ID}.bam

# Calculate statistics
samtools flagstat $OUTPUT_DIR/${SAMPLE_ID}.bam > $OUTPUT_DIR/${SAMPLE_ID}.flagstat

# End timing
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Job completed in $elapsed_time seconds"
```

