# Simple Command-Line Tutorial for Conda-Based Tools

A **quick guide** to using **Conda**, **Anaconda**, **Miniconda**, **Miniforge**, **Mamba**, and **Micromamba**.

Here is a simple command-line tutorial overview for the most common Conda-based tools: Conda, Anaconda, Miniconda, Miniforge, Mamba, and Micromamba.

---

## 1. Conda (General Overview)

`conda` is a **package and environment manager**, forming the backbone of Anaconda, Miniconda, etc. It allows you to:  
 - Create and manage **isolated environments** for different projects  
 - Install specific versions of tools/packages  

### Installing conda
There are three ways/options for installing conda:

**Option 1:** Install Conda via Anaconda (Full Package)   
**Option 2:** Install Conda via Miniconda (Minimal Package)  
**Option 3:** Install Conda via Miniforge (Community Version)  

See sections below for a more detailed explanation.  
Once Conda is installed, you can create environments, install packages, and manage dependencies easily!  

## 2. Anaconda

- **Anaconda** = Conda + 100s of preinstalled data science packages (~3–5 GB).
- Ideal for **full-featured local data science** setups.

### Install:

Download from: [Anaconda official website](https://www.anaconda.com)

### Usage:

Once installed, you use `conda` commands as shown above—no difference.

---

## 3. Miniconda

- A **minimal** version of Anaconda with just `conda`, Python, and essential packages.
- Ideal for creating custom environments from scratch.

### Install:

Download from: [Miniconda official site](https://docs.conda.io/en/latest/miniconda.html)

### Usage:

Same as `conda` commands above.

---

## 4. Miniforge

- Like Miniconda but **built using open-source conda-forge packages**.
- No proprietary packages.
- Recommended if you use **only `conda-forge` channel**.

### Install:

[Miniforge GitHub](https://github.com/conda-forge/miniforge)

```bash
# Example: Create env using conda-forge
conda create -n myenv -c conda-forge python=3.11
```

---

### Basic Conda Commands

```bash
# Check conda version
conda --version

# Update conda
conda update conda

# Create a new environment
conda create -n myenv python=3.11

# Activate the environment
conda activate myenv

# Deactivate current environment
conda deactivate

# Install a package
conda install numpy

# List installed packages
conda list

# List all environments
conda env list

# Remove environment
conda remove -n myenv --all
```

---

## 5. Mamba

- A **faster drop-in replacement** for Conda (written in C++).
- Speeds up solving and installing environments significantly.

### Install into Base Environment:

```bash
conda install mamba -n base -c conda-forge
```

### Usage (same as Conda but replace `conda` with `mamba`):

```bash
mamba create -n fastenv python=3.11
mamba install numpy pandas
```

---

## 6. Micromamba

- A **lightweight, standalone version of Mamba**. No Python required.
- Ideal for **Docker containers** or minimal environments.

### Install (Linux/macOS example):

```bash
curl micro.mamba.pm/install.sh | bash
# Add to shell: eval "$(~/micromamba/bin/micromamba shell hook -s bash)"
```

### Usage:

```bash
micromamba create -n myenv python=3.11
micromamba activate myenv
micromamba install numpy
```

---

## Summary Table

| Tool       | Size         | GUI | Use Case                              |
| ---------- | ------------ | --- | ------------------------------------- |
| **Anaconda**   | Huge (~5GB) | ✅ Yes | Full data science suite               |
| **Miniconda**  | Small        | ❌ No | Minimal + customizable Conda          |
| **Miniforge**  | Small        | ❌ No | Open-source alternative (`conda-forge`) |
| **Conda**      | N/A          | N/A | The tool powering the others          |
| **Mamba**      | Add-on       | ❌ No | Faster Conda                          |
| **Micromamba** | Minimal      | ❌ No | Lightweight, ideal for containers     |

---

## Exercises: Now lets embark on creating and installing bioinformatics softwares in an environment.  

**Problem Statement:** *Using mamba create an environment called 'acdc_amr', installed with python verion 3.11 and activate it, install fastp, bowtie, spades and megahit all from bioconda channel. Then export to a yaml file. The ymal file will ensure* **portability** and **reproducibility**

**Solution:**
Creating a **reproducible environment** with `mamba` ensures consistency across different systems, allowing easy replication of analyses. Here’s a **step-by-step guide** to setting up and exporting a reproducible Conda environment using `mamba`:

---

## **Step 1: Install Mamba**  

If you don’t have `mamba`, install it into your base Conda environment:
```bash
conda install mamba -n base -c conda-forge
```
  > `mamba` is a **drop-in replacement** for `conda` but much faster in resolving dependencies.

---

## **Step 2: Create an Environment Called `acdc_amr` with Python 3.11**  

Use the `mamba create` command to set up the environment:

```bash
mamba create -n acdc_amr python=3.11
```
 >  `-n acdc_amr`: Names the environment `acdc_amr`    
 >  `python=3.11`: Specifies Python version **3.11**  

---

## **Step 3: Activate the Environment**

After creation, activate the new environment:

```bash
mamba activate acdc_amr
```

 > This switches from the **base environment** to `acdc_amr`.

---

## **Step 4: Install Required Packages from `bioconda`**  

Now, install **fastp**, **bowtie**, **SPAdes**, and **megahit** using `mamba`, making sure they come from the `bioconda` channel:  

```bash
mamba install -c bioconda fastp bowtie spades megahit
```

 > `-c bioconda` → Ensures all packages are fetched from **bioconda**  
 > This installs sequencing tools critical for bioinformatics analysis.

---

## **Step 5: Export Environment to a YAML File**  

To ensure the environment can be **reproduced** exactly, export its configuration:

```bash
mamba env export -n acdc_amr > acdc_amr_env.yaml
```

**Creates a file** (`acdc_amr_env.yaml`) containing:
- **Python version**
- **All installed packages**
- **Channels used (`bioconda`)**

---

## **Step 6: Recreate the Environment from YAML File**

Anyone can **replicate the environment** using:

```bash
mamba env create -f acdc_amr_env.yaml
```

 Ensures that another user gets **exactly the same environment**!

---

### ** Why is this Important for Reproducibility?**
✔ **Version Control** → Fixes Python and package versions  
✔ **Channel Control** → Ensures all dependencies come from `bioconda`  
✔ **Portable** → YAML file allows easy recreation across machines  


## Next Steps

Would you like to know alternative package management system for a HPC?  
See the next part on HPC module management system:

[Introduction to Module System](https://github.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/blob/main/Intro_Module_System.md)

Happy coding! 🎯

