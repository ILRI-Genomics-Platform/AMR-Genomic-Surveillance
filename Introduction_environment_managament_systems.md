# Simple Command-Line Tutorial for Conda-Based Tools

A **quick guide** to using **Conda**, **Anaconda**, **Miniconda**, **Miniforge**, **Mamba**, and **Micromamba**.

Here is a simple command-line tutorial overview for the most common Conda-based tools: Conda, Anaconda, Miniconda, Miniforge, Mamba, and Micromamba. markdown document

---

## 1. Conda (General Overview)

`conda` is a **package and environment manager**, forming the backbone of Anaconda, Miniconda, etc.

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

## Next Steps

Would you like to know alternative package management system for a HPC?
See the next part on HPC module management system:

[Introduction to Module System](https://github.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/blob/main/Intro_Module_System)

Happy coding! 🎯

