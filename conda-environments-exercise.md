# Conda environment exercise
To practice and learn more about `conda` and its environments, work on the following exercise by following the instructions below.

1. Having installed `conda` through Miniforge, check the version of `conda` that you have.
  - Additionally, check the version of `mamba`
2. Create a `conda` environment called: fastqc-env
3. Activate the `conda` environment that you created in step 2 above.
4. Once activated, install the bioinformatics tool `fastqc` in the environment using `conda`.
5. Check the version of `fastqc` that you have installed.
6. Create a directory called `conda_exercise` in your home folder and change directory to the created folder.
8. Use the Linux tool `wget` to download the file from this url: https://hpc.ilri.cgiar.org/~kmwangi/SRR25008772_1.fastq.gz
9. Run the following command `fastqc SRR25008772_1.fastq.gz`. Can you explain what the command does?
10.Download the following `yml` file that contains conda instruction on installing certain bioinformatics tools: https://hpc.ilri.cgiar.org/~kmwangi/qc-tools-environment.yml
  - Which tools will be installed if that `yml` file is used as an argument for conda?
  - Into which environment will be the tools be installed?
  - From which conda channels are we obtaining the tools?
11.Use the downloaded yml file in step 11 above to install the tools.
12.On your terminal, list all the conda environments that you have. Can you identify the active environment?

