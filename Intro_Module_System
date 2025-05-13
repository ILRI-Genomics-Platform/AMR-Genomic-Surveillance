Here is a **simple command-line tutorial** for using the **module system**—a common tool for managing software on **HPC (High-Performance Computing)** systems.

---

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

## Advanced Tip: Persistent Modules

To **auto-load modules** on login, add `module load` lines to `~/.bashrc`, `~/.bash_profile`, or an HPC-specific profile script:

```bash
echo "module load python/3.10" >> ~/.bashrc
```

---

