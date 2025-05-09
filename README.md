# Bioinformatics Analysis for Antimicrobial Resistance Genomic Surveillance  

---  

###### **_Trainers_**: [John Juma](https://github.com/ajodeh-juma), [Kennedy Mwangi](https://github.com/wanjauk), [Ouso Daniel](https://github.com/ousodaniel) & [Gilbert Kibet](https://github.com/kibet-gilbert)

---

- **Table of Contents**
  - [Introduction](#introduction)
  - [Scope of the tutorial](#scope-of-the-tutorial)
  - [pre-requisites](#pre-requisites)
  - [Set up](#set-up)
  - [Analysis preparations](#analysis-preprations)
    - [Logging into the HPC](#logging-into-the-hpc)


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

From here, we move to specific tutorial for analysis.
