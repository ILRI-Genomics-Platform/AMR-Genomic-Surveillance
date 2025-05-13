# Conda and bioinformatics software environment management
In this tutorial, you will learn how to create a conda environment and install tools in the environment.
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
conda activate fastqc-env
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

