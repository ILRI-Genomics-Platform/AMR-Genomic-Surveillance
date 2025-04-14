# Bioinformatics Analysis for Antimicrobial Resistance Genomic Surveillance  

---  

###### **_Trainers_**: [John Juma](https://github.com/ajodeh-juma), [Kennedy Mwangi](https://github.com/wanjauk) & [Gilbert Kibet](https://github.com/kibet-gilbert)

---

- **Table of Contents**
  - [Introduction](#introduction)
  - [Scope of the tutorial](#scope-of-the-tutorial)
  - [Background](#background)
  - [pre-requisites](#pre-requisites)
  - [Set up](#set-up)
  - [Analysis preparations](#analysis-preprations)
    - [Logging into the HPC](#logging-into-the-hpc)
    - [Project organisation](#project-organisation)
    - [Data retrieval](#data-retrieval)
  - [Bioinformatics analysis](#bioinformatics-analysis)
    - [ONT](#ont)
    - [Illumina](#illumina)
      - [Step 1: Loading modules](#step-1-loading-modules)
      - [Step 2: Data quality assessment](#step-2-data-quality-assessment)
      - [Step 3: Quality and adapater filtering](#step-3-quality-and-adapater-filtering)
      - [Step 4: Reads alignment to reference genome](#step-4-reads-alignment-to-reference-genome)
      - [Step 5: Generate genome consensus](#step-5-generate-genome-consensus)
      - [Step 6: De novo genome assembly](#step-6-de-novo-genome-assembly)
      - [Step 7: Quality assessment of the assembled contigs](#step-7-quality-assessment-of-the-assembled-contigs)
      - [Step 8: Annotatio of the assembled contigs](#step-8-annotatio-of-the-assembled-contigs)
      - [Step 9: Identification of AMR genes](#step-9-identification-of-amr-genes)
      - [Step 10: Identification of the virulence factors](#step-10-identification-of-the-virulence-factors)
  - [References](#references)

## Introduction  
*Escherichia coli* are ubiquitous bacteria found in the environment including the gut of humans and other animals and consists of numerous types and strains. Most types of *E. coli* are normal inhabitants of the gut of animals and do not cause diseases. However, some of the strains of *E. coli* have acquired genes that enable them to cause disease. These strains are commonly associated with food poisoning leading to diarrhoea and are referred to as diarrheagenic *E. coli* (DEC). Transmission occurs primarily through contaminated food but can also occur via person-to-person transmission, animal contact and water. For example, Shiga toxin-producing *E. coli* (STEC) serotype O157:H7 causes bloody diarrhoea and has previosuly been responsible for outbreaks worldwide.   

*Klebsiella pneumoniae* -  

## Scope of the tutorial  
In this workshop we will tackle, hands-on, the basic principles employed to generate consensus genome sequences of *E. coli* and *K. pneumoniae* and identify the serotypes involved in outbreaks, anti-microbial resistance genes (AMRs) and possible virulence factors they may have. 

## Background  


## Pre-requisites  
This module will come after the Introductory Linux module and therefore assumes familiarity with basic Linux command-line use. It also assumes you have an account and are operating in the ILRI computing cluster from a local Linux environment. 

## Set up  
You will use your personal computers computers to log into the ILRI HPC cluster, which operates on a Linux-operating system. Since we will be working from the remote servers, you will not need special setup for your personal laptops. However, you will need to install a program that enables you to log into the HPC.

>**Note**

>Once inside the `hpc`, all instances of ```$USER``` will be equivalent to the hpc username that you were assigned. Your username, by default, is stored in a variable called `USER`. By using it, you will not have to type-in your username, rather, your shell will automatically pick your username which is the value stored in the `USER` variable. The `$` (dollar) character-prefix to a variable name is used to call the value of that variable.



## Analysis preprations

### Logging into the HPC  
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

### Project organisation  
We will start by setting up the project directory structure and then conduct the analysis stepwise. To setup a well-structured project directory we need to create some directories to store our data and scripts. We will be conducting our a anlysis from a directory in the `scratch` space of the HPC.  

1. *Create a directory using your username in the scratch:*
> **Note:** In this command we use the UNIX environment variable `$USER` which by default was created during logging into the HPC to store your `<user_name>` i.e (`Bio4InfoXX`). You can view its value using the command `echo $USER`.  
```
mkdir -p /var/scratch/$USER
cd /var/scratch/$USER
```
2. *Create project directories:*
> **Note:** We create a project directory `viralMetagen` to store all that pertains to this tutorial/project. Within `viralMetagen` we created `data` and subdirectories to store our input data and results from different analysis steps. We create `scripts` directory to store scripts/code that we genenrate or need in the analysis.
```
mkdir -p amr-surveillance/{data,scripts}
cd amr-surveillance/
mkdir -p ./data/{database,fastq,fastqc,fastp,spades,quast,bowtie,samtools}
```


### Data retrieval  

## Bioinformatics analysis  

### ONT 


### Illumina  

### Step 1: Loading modules  
### Step 2: Data quality assessment 
### Step 3: Quality and adapater filtering  
### Step 4: Reads alignment to reference genome    
### Step 5: Generate genome consensus  
### Step 6: De novo genome assembly  
### Step 7: Quality assessment of the assembled contigs  
### Step 8: Annotatio of the assembled contigs   
### Step 9: Identification of AMR genes  
### Step 10: Identification of the virulence factors  

## References  
