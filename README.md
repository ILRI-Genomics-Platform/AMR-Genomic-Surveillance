# Bioinformatics Analysis for Antimicrobial Resistance Genomic Surveillance  

---  

###### **_Trainers_**: [John Juma](https://github.com/ajodeh-juma), [Kennedy Mwangi](https://github.com/wanjauk), [Ouso Daniel](https://github.com/ousodaniel) & [Gilbert Kibet](https://github.com/kibet-gilbert)

---

- **Table of Contents**
  - [Introduction](#introduction)
  - [Scope of the tutorial](#scope-of-the-tutorial)
  - [Pre-requisites](#pre-requisites)
  - [Set up](#set-up)
  - [Analysis preparations](#analysis-preprations)
    - [Required tools](#required-tools)
    - [Logging into the HPC](#logging-into-the-hpc)
  - [Project organisation](#project-organisation)
  - [Software installation (local computer)](#software-installation-local-computer)
    - [Creating a conda environment](#creating-a-conda-environment)
    - [Creating a Python environment](#creating-a-python-environment)
  - [Analysis tutorials](#analysis-tutorials)


## Introduction  
*Escherichia coli* are ubiquitous bacteria found in the environment including the gut of humans and other animals and consists of numerous types and strains. Most types of *E. coli* are normal inhabitants of the gut of animals and do not cause diseases. However, some of the strains of *E. coli* have acquired genes that enable them to cause disease. These strains are commonly associated with food poisoning leading to diarrhoea and are referred to as diarrheagenic *E. coli* (DEC). Transmission occurs primarily through contaminated food but can also occur via person-to-person transmission, animal contact and water. For example, Shiga toxin-producing *E. coli* (STEC) serotype O157:H7 causes bloody diarrhoea and has previosuly been responsible for outbreaks worldwide.   

*Klebsiella pneumoniae* is a Gram-negative, non-motile, encapsulated, lactose-fermenting and facultative anaerobic, rod-shaped bacterium.  It can be found in a wide range of environments from soil to plants, water or human bodies where it is a normal flora of mucosal membranes of the mouth, skin, and intestines. However, once it enters the human body or animal lungs, it becomes pathogenic; highly virulent and resistant to antibiotics. It is the most common cause of hospital-acquired pneumonia and a very important cause of nosocomial (hospital-acquired) bacterial infections. It has also been associated with patients with alcohol use disorder or diabetes mellitus but also a major threat to newborns, elderly and immunocompromised patients.

## Scope of the tutorial  
In this workshop we will tackle, hands-on, the basic principles employed to generate consensus genome sequences of *E. coli* and *K. pneumoniae* and identify the serotypes involved in outbreaks, anti-microbial resistance genes (AMRs) and possible virulence factors they may have. 

## Pre-requisites  
This module will come after the Introductory Linux module and therefore assumes familiarity with basic Linux command-line use. It also assumes you have an account and are operating in the ILRI computing cluster from a local Linux environment. 

## Set up  
You will use your personal computers computers to log into the ILRI HPC cluster, which operates on a Linux-operating system. Since we will be working from the remote servers, you will not need special setup for your personal laptops. However, you will need to install a program that enables you to log into the HPC.

## Analysis preprations

### Required tools
If you are using a Windowas OS computer, you will need to install the follwing tools to enable you access the remote HPC server. Additional tools will help
in downloading data from HPC to your local computer.

Please follow the instructions here: [List of required tools](https://github.com/ILRI-Genomics-Platform/trainings-required-software)

### Logging into the HPC
On your computer open a terminal. On Windows use Windows Subsystem for Linux (WSL) or MobaXterm. If you are using Linux, click on the Terminal icon to open the terminal.  
To sign in to HPC using the following command, you will have been assigned a ***username*** that looks like `Bio4InfoXX` and a ***password***.
1. Replace `<user_name>` in the command with the provided username and execute (click enter). 
2. Enter the password and execute. ***Note:*** The password will be typed in the background but will not be visible to you as you type.
```
ssh <user_name>@hpc.ilri.cgiar.org
```
There are two nodes to choose from: `compute05`  and `compute06`. If your username (`Bio4InfoXX`) ends with an ***Odd Number*** (1,3,5,7,9) use `compute05` and if it ends with n ***even number*** (2,4,6,8,0) use `compute06`. Now let us secure a four of CPUs in one of the HPC nodes.  
>Compute05
```
interactive -w compute05 -c 2 -J amr-surveillance -p batch
```
>Compute06
```
interactive -w compute06 -c 2 -J amr-surveillance -p batch
```

# Project organisation  
We will start by setting up the project directory structure and then conduct the analysis stepwise. To setup a well-structured project directory we need to create some directories to store our data and scripts. We will be conducting our a anlysis from a directory in the `scratch` space of the HPC.  

1. *Create a directory using your username in the scratch:*
>**Note**

>Once inside the `hpc`, all instances of ```$USER``` will be equivalent to the hpc username that you were assigned. Your username, by default, is stored in a variable called `USER`. By using it, you will not have to type-in your username, rather, your shell will automatically pick your username which is the value stored in the `USER` variable. The `$` (dollar) character-prefix to a variable name is used to call the value of that variable.

```
mkdir -p /var/scratch/$USER
cd /var/scratch/$USER
```
2. *Create project directories:*
> **Note:** 

> We create a project directory `ACDC_AMR2025` to store all that pertains to this tutorial/project. Within `ACDC_AMR2025` we will have `data` and subdirectories to store our input data and `results` from different analysis steps. We will also have `scripts` directory to store scripts/code that we genenrate or need in the analysis.

```
mkdir -p ACDC_AMR2025
cd ACDC_AMR2025
```

From here, we move to specific tutorial for analysis.

## Software installation (local computer)
Please follow the following instructions to install the tools that we will use on local computer. Note that these tools have already been installed on the remote computer (the HPC) by the system admin.

Open a terminal in using WSL in Windows. If you have a Linux OS, open a terminal by clicking on the terminal icon. 

It is good practice to organize your work in project directories so as to have all files relating to a project in one place. Therefore, let's create a directory
on our local machine. We will create a directory similar to the one we created on the HPC above. Run the following command:

```
mkdir -p ACDC_AMR2025
cd ACDC_AMR2025
```

Download the following `yml` file to your computer.
```
wget https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/environment.yml
wget https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/snippy-environment.yml
wget https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/artic-mpox-environment.yml
wget https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/artic-requirements.txt
```

## Analysis tutorials
Please find here linked the tutorials for Introduction to Linux and for respective organisms. Right click on the link and open it in a new tab on your browser.

1. [Introduction to Linux](linux_hpc_intro.md)
2. Environment Management Systems
   - [Setting up `Conda`](Setting_up_conda.md)
   - [Environment Manament with Conda](Introduction_environment_managament_systems.md)
   - [Introduction to HPC Module System](Intro_Module_System.md)
   - [Running Jobs in the HPC using SLURM](Introduction_to_High_Performance_Computing_with_SLURM.md)
3. [*Klebsiella pneumoniae*](bacterial-amr-analysis-for-ont.md) - Long Reads (ONT)
   - [Introduction to *K. pneumoniae*](Intro_Klebsiella_pneumoniae.md)
   - [*K. pneumoniae* Bacterial Whole Genome Assembly and AMR analysis](bacterial-amr-analysis-for-ont.md)
4. [*Escherichia coli*](bacterial-amr-analysis-for-illumina.md) - Short Reads (Illumina)
   - [Introduction to *E. coli*](Intro_Escherichia_coli.md)
   - [*E. coli* Bacterial Whole Genome Assembly and AMR analysis](bacterial-amr-analysis-for-illumina.md)
5. [*Monkeypox virus*](mpox-genome-analysis-illumina.md)
   - [Introduction to *Monkeypox Virus* (MPXV)](Intro_Monkeypox_Virus.md)
   - [Monkeypox Virus Whole Genome Sequene analysis](mpox-genome-analysis-illumina.md) - Short Reads (Illumina)
6. [Bioinformatics Data Management](Intro_to_bioinformatics_data_Management.md)


> [Bacterial Genome and AMR Analysis Terms](bacterial_genome_terms.md)
