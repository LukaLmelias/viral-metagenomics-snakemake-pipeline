#!/bin/bash

### DAILY TRACK OF COMMANDS

### 15/05/2023: on my birthday :)

#ls

#setting alias to luster:
#nano ~/.bashrc

#alias lustre_luka="cd /lustre/shared/wfsr-mcfa/projects/internships/luka"

#source ~/.bashrc


#######################  INSTALLING TOOLS

#snakemake
# conda install -c bioconda snakemake # this never worked ):


# mamba create -c conda-forge -c bioconda -n snakemake snakemake # this worked


#installing fastqc
#conda install -c bioconda fastqc

#installing trimmomatic
#conda install -c bioconda trimmomatic



# installing multiqc
#conda install -c bioconda multiqc  #failed ):

#some tools only work better with pip
# so install pip within conda


#installing megahit
#conda install -c bioconda megahit  


#install metaSPAdes; actually it is part of SPAdes
#conda install -c bioconda spades

#installing kraken2
#conda install -c bioconda kraken2 # this conflicts could not be resolved


#install meta
#conda install -c bioconda metaquant



#the conflict have become too many to resolve 
# i will remove the original conda environmwent and set up each tool it own environment
# hopefully this will avoid more conflicts

#conda remove -n snakemake --all

# NOW INSTALL EVERYTHING ON ITS OWN

# here is the plan:
    #conda create the tool env
    #conda install the tool
    #conda activate the tool env
    #test that the tool properly installed



# snakemake 
#mamba create -c conda-forge -c bioconda -n snakemake snakemake


# i decided to write a tiny script to install them on by one
# tools[ fastqc, multiqc, trimmomatic,
#megahit, spades, metaquant, kraken2]

# echo hi, tool name?:
# read tool
# mamba create -c conda-forge -c bioconda -n $tool $tool 



#next is to add them on the environment path
# it aint working ):

#####################################
#####################################

### 16/06/2023

#####################################
#####################################

# running conda envs in snakemake

# export the conda environments

# read env_name
# conda activate  $env_name

# conda env export > $env_name.yml

# REMEMBER TO USE CROSSPLATFORM PINNEDE

# alternative is to use docker




#Some deleted workflows



# rule hello:
#     input:
#         R1 = expand("{dir}{forward}", dir=raw_reads_path,forward=FORWARD),
#         #R2 = expand("{dir}{reverse}", dir=raw_reads_path,reverse=REVERSE),
#     output:
#         out1= "{sample}_1.txt",
#         #out1 = expand("{sample}_1.txt", sample=SAMPLES),

#     shell:
#         """
#         head -n 1  {input.} > {output.out1}
        

#         """

# rule fastp:
#     input: 
#         R1 = expand("{dir}{forward}", dir=raw_reads_path,forward=FORWARD),
#         R2 = expand("{dir}{reverse}", dir=raw_reads_path,reverse=REVERSE),
#     output:
#         "output.txt"
    
#     conda:
#        "/lustre/shared/wfsr-mcfa/projects/internships/luka/tools/yamls/fastp.yml"
#     shell:
#        "fastp -? > {output}"

