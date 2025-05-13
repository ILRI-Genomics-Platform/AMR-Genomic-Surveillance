# **A Simple Command-Line Tutorial** for using the *module system*

**module system** is a common tool for managing software on **HPC (High-Performance Computing)** systems.


## What is the Module System?

The **Environment Modules System** lets you easily **load, unload, and switch** between software versions without manually editing paths.

Used on many HPC clusters to manage software.

---

## Common Module Commands

### 1. **List Available Modules**

```bash
module avail
```

Shows all software packages and versions available to load.

---

### 2. **Load a Module**

```bash
module load <name>/<version>
```

Example:

```bash
module load fastp/0.24.1
module load python/3.10
```
Alternatively:  
```
module load fastp/0.24.1 python/3.10
```

You can omit the version if the default is fine:

```bash
module load python
```

---

### 3. **Unload a Module**

```bash
module unload <name>
```

Example:

```bash
module unload python
```

---

### 4. **Swap Modules**

```bash
module swap <old> <new>
```

Example:

```bash
module swap fastp/0.24.1 fastp/0.22.0
```

---

### 5. **List Loaded Modules**

```bash
module list
```

---

### 6. **Unload All Modules**

```bash
module purge
```

---

### 7. **See Module Help**

```bash
module help <name>
```

Example:

```bash
module help python
```

---

## Example Workflow

```bash
# Check available Python versions
module avail python

# Load Python 3.11
module load python/3.11.1

# Verify it's active
python --version

# Load a Reads qc
module load fastp/0.24.1

# Swap to another fastq reads qc
module swap fastp/0.24.1 fastp/0.22.0

# Clean up environment
module purge
```

---
## Pros and X Cons of the Module System:  


|✅ Pro                                        | ❌ Cons                                 |
| -------------------------------------------- | --------------------------------------- |
| Easy switching between software versions     | Can get messy without structure         |
| Per-user customization                       | Manual setup for custom-built software  |
| Reduces path conflicts                       | Not a full package manager (no install) |
| Works well with job schedulers (e.g., SLURM) | Requires sysadmin effort to set up      |
| Clean environment separation                 | Harder to track dependencies            |

---
## Advanced Tip: Persistent Modules

To **auto-load modules** on login, add `module load` lines to `~/.bashrc`, `~/.bash_profile`, or an HPC-specific profile script:

```bash
echo "module load python/3.10" >> ~/.bashrc
```

---
## Execises: We will now try to load differnt softwares we may need in the AMR analysis and resolve conflicts.

**Problem Statement:** You are preparing to execute some analysis steps for bacterial AMR based on ONT data. Here are some of the packages that will be needed: `samtools/1.9`, `racon/1.5.0`, `unicycler/0.4.7`, `prodigal/2.6.3`, `fastp/0.22.0`, `porechop/0.2.4`, `minimap2/2.13`, `spades/3.13.0`, `mlst/2.23.0`, `infernal/1.1.2`, `fastqc/0.11.9`, `any2fasta/0.4.2`, `medaka/0.8.2`, `velvet/1.2.10`, `hmmer/3.3`, `prokka/1.14.6`, `lighter/1.1.2`,  `flye/2.4.2`, `megahit/1.2.9`, `bowtie2/2.3.4.1`, `bedtools/2.29.0`,  `flash/1.2.11`, `miniasm/0.3`, `nanoplot/1.42.0`, `nanoq/0.10.0` and `porechop`. Please load them using `module load` command, identify dependancies and conflicts where there are any.

### **Step 1:**  Attempt to load all packages at onces.  

```bash
module loadd samtools/1.9 racon/1.5.0 unicycler/0.4.7 prodigal/2.6.3 fastp/0.22.0 porechop/0.2.4 minimap2/2.13 spades/3.13.0 mlst/2.23.0 infernal/1.1.2 fastqc/0.11.9 any2fasta/0.4.2 medaka/0.8.2 velvet/1.2.10  hmmer/3.3 prokka/1.14.6 lighter/1.1.2 flye/2.4.2 megahit/1.2.9 bowtie2/2.3.4.1 bedtools/2.29.0 flash/1.2.11 miniasm/0.3  nanoplot/1.42.0 nanoq/0.10.0 porechop/0.3.2pre
```

 - Which dependancies did you see?  
 - What are the conflicts you noticed?  
 - How best do you think you can resolve these conflicts?  
 - Why are there conflicts with `porechop`?  

### **Step 2:**  List modules that have been loaded.  

```bash
module list
```

 - Are there extra modules that have been loaded? and Why?  
 - Can you swap `porechop` to latest version?  

### **Step 3:**  Lets manage our modules even further...  

 - `module avail` all versions of `samtools`. Which version is loaded and which version is the latest?  
 - Can you `module unload` current version and `module load` latest version?  
 - `module swap` the current version to another version of `samtools` 
 > Now `module purge` all modules and check if you succeeded.   

### **Step 4:**  Lets now discuss the pros and cons of **Module System** vs **conda** System.  

 - Are there any benefits of module system vs conda system in a PC vs a HPC?  
 - How easy is each to set-up/install?  
 - What would you wish to work with in your PC vs a shared cluster?  

