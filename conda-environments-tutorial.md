# Conda and bioinformatics software environment management

## Installing Miniforge
Please follow the instructions below which are obtained from miniforge [GitHub repository](https://github.com/conda-forge/miniforge?tab=readme-ov-file#unix-like-platforms-macos-linux--wsl)

### Linux-based distributions
Download the installation bash script.
```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
```

Run the script to install Miniforge.
```
bash Miniforge3-$(uname)-$(uname -m).sh
```
To add conda to your `.bashrc` for convenient initialization any time you log in to the HPC, please run the following code:
```
eval "$(/home/${USER}/miniforge3/bin/conda shell.bash hook)"

conda init 

```

### MacOS
If you are using MacOS, install Miniforge by first installing Homebrew:
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Once already installed, install Miniforge as follows:
```
brew install miniforge
```

Finally run this command:
```
conda init zsh
```
Restart your terminal.

## Check your conda installation

We will now learn how to create a conda environment and install tools in the environment.
This tutorial assumes that you have already installed one of the minimal installers of `conda` such as `Miniconda` or `Miniforge`.

We can check that `conda` is installed by printing out its help page or by printing out its version
```
conda --version
```

```
conda --help
```

## Creating an environment
First, list environments to see what you get. Run the following code:

```
conda info --envs
```

Now, create an environment called `fastp-env` as follows:

```
conda create -n fastp-env
```

List the environments again to see what was added to the list of environments.

Activate the environment that you created as follows:

```
conda activate fastp-env
```

Did you see any change on your terminal?

Now change the environment back to your default one as follows:

```
conda deactivate
```

## Installing packages to an environment

First activate the environment that you created. You can check the tool that you want to install on the anaconda webpage to find out more information 
about it. For example, here is the webpage for `fastp` showing information about its version and how to download it: https://anaconda.org/bioconda/fastp 
You can search any tool in the search bar.

Next, install the `fastp` tool as follows:
```
conda install bioconda::fastp
```

Here we see that we obtained the `fastp` tools from the `bioconda` channel. 

We can also specicfy the version of `fastp` that we want to downlaod as follows:

```
conda install bioconda::fastp=0.24.1
```
Another way to install the tools is by listing them in a `yaml` file and then passing that file as an argument to `conda`.
Download the following `yml` file and use it to create and install `fastp` in another environment.

```
wget https://hpc.ilri.cgiar.org/~kmwangi/environment.yml
```

Install the `fastp` tool using `conda` as follows:
```
conda env create -f environment.yml
```

