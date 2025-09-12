# Reproducing the analysis environments with *mamba*

You’re about to set up your bioinformatics analyis tools with **Mamba**, a faster, smarter alternative to Conda.
Before we do let us download the enviroment recipe files i.e the YAML files.

---

## What Is a Conda/Mamba Environment YAML File?

A **YAML file** for Conda or Mamba is a **blueprint** for creating a reproducible software environment. It defines:

- **Channels**: Where to fetch packages (e.g., `bioconda`, `conda-forge`)
- **Dependencies**: The exact tools and versions you need
- **Optional pip packages**: For Python tools not available via Conda

This is especially useful in bioinformatics, where tool conflicts are common — like in your case with ONT and Illumina pipelines.

---

### Example Breakdown

```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - bioconda::fastqc=0.11.9
  - bioconda::prokka=1.14.6
  - pip
  - pip:
    - nanoplot==1.42.0
```

This YAML:
- Pulls packages from `bioconda` and `conda-forge`
- Installs `fastqc`, `prokka`, and `nanoplot`
- Ensures version control for reproducibility

---

## Download YAML Files from GitHub

The YAML files for these tutorials are hosted in GitHub in our course repo. Here’s how to download them:

### Option 1: Using `wget`

```bash
wget https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/envs/env_ont-Ill_klebStep1-6.yaml
wget https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/envs/env_ont_klebStep7-amr.yaml
wget https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/envs/env_ont_klebStep8-10.yaml
wget https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/envs/mpox-environment.yaml
```

- `env_ont-Ill_klebStep1-6.yaml` - Has packages for Step 3 to Step 6 for both ONT and Illumina data.
- `env_ont_klebStep7-amr.yaml` - Has packages for Step 7 AMR prediction using ResFinder, RGI and AMRFinderPlus.
- `env_ont_klebStep8-10.yaml` - Covers Step 8 and 9 packages.
- `mpox-environment.yaml` - Has all packages needed by the Mpox workflow. That also includes packages installed in `Step 2`, i.e packages listed in `./utils/artic-requirements.txt` so skip `Step 2`.


### Option 2: Using `curl`

```bash
curl -O https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/envs/env_ont-Ill_klebStep1-6.yaml
```

---

## Why Separate Environments?

Since ONT and Illumina tools for different Steps can conflict (e.g., different versions of `samtools`, `minimap2`, or `shovill`), using **three isolated environments** lets you:

- Switch between pipelines without breaking dependencies
- Run analyses in parallel
- Keep your environments lean and purpose-built

---

Now let us proceed with the environment setup by starting with checking if `conda` and `mamba` are installed.

---

## Step 1: Check if Conda is Installed

Open your terminal and run:

```bash
conda --version
```

If you see something like `conda 24.3.0`, you're good. If not, you’ll need to install Conda first — typically via [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Miniforge3](https://github.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/blob/main/Setting_up_conda.md) for a lightweight setup.

---

## Step 2: Check if Mamba is Installed

Run:

```bash
mamba --version
```

If you get a version number, Mamba is already installed. If you see `command not found`, it’s time to bring in the upgrade.

---

## Why Mamba > Conda

Mamba is a **drop-in replacement** for Conda, but it’s built in C++ and optimized for speed. Here’s why it’s better:

| Feature         | Conda            | Mamba             |
|----------------|------------------|-------------------|
| Speed          | Slower dependency resolution | ⚡ Blazing fast |
| Parallelism    | Single-threaded  | Multi-threaded    |
| Reliability    | Occasional solver failures | More robust solver |
| Syntax         | Same as Conda    | Same as Conda     |

If you’ve ever waited minutes for Conda to resolve dependencies, Mamba will feel like a breath of fresh air.

---

## Step 3: Install Mamba

You can install Mamba into your base Conda environment:

```bash
conda install mamba -n base -c conda-forge
```

This installs Mamba from the `conda-forge` channel and makes it available globally.

---

## Step 4: Create Environment Using Mamba

Assuming you have a YAML file like `env_ont-Ill_klebStep1-6.yaml`, run:

```bash
mamba env create -f ./env_ont-Ill_klebStep1-6.yaml --name bactWGSgenome
```

---

## Pro Tip: Validate Before You Build

To preview what Mamba will install:

```bash
mamba env create -f ./env_ont-Ill_klebStep1-6.yaml --name bactWGSgenome --dry-run
```

This helps catch issues before committing.

---

## Step 5: Activate, Use and Deactivate Environment

Then activate it:

```bash
mamba activate bactWGSgenome
```

And deactivate when done:

```bash
mamba deactivate
```

---

For the other two environments repeat **Step 4 and 5** for each YAML file.
Example:

```bash
# bactWGSamr environment
mamba env create -f env_ont_klebStep7-amr.yaml --name bactWGSamr --dry-run
mamba env create -f env_ont_klebStep7-amr.yaml --name bactWGSamr
mamba deactivate

# bactWGSvc environment
mamba env create -f env_ont_klebStep8-10.yaml --name bactWGSvc --dry-run
mamba env create -f env_ont_klebStep8-10.yaml --name bactWGSvc
mamba activate bactWGSvc
mamba deactivate

# mpoxWGSvc environment
mamba env create -f mpox-environment.yaml --name mpoxWGSvc --dry-run
mamba env create -f mpox-environment.yaml --name mpoxWGSvc
mamba activate mpoxWGSvc
mamba deactivate
```

Beyond the virtual environments, there are `python` and `R` scripts used in different steps of the pipeline. For example the mpox pipeline uses these python scripts in different steps:  
 - `align_trim.py` - `Step 7`,   
 - `process_gvcf.py` - `Step 8` and  
 - `fetch-genomes.py` - `Step 10`  

Some of these scripts can be found in the [**scripts**](../scripts) directory of these GitHub repo. To access them;  
 - Go to the repo's [scripts](../scripts) ->   
 - Click on the script e.g [fetch-genomes.py](https://github.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/blob/main/scripts/fetch-genomes.py) -> 
 - Click the `Raw` button, copy the URL (link) and download using `wget <URL>` as shown **OR** 
   - `Copy raw file` and paste it to a file with same name **OR** 
   - click `Download raw file` and move it to the right folder.  

