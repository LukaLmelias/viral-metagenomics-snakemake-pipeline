#!/bin/bash


# """
# This script installs snakemake 
# creates snakemake conda/mamba environment environment
# It assumes conda/minin-conda is installed

# Before running this script, edit the config file paths (have a separate script to create configs)

# """

# check to install mamba (check if its installed)
#from:https://stackoverflow.com/questions/592620/how-can-i-check-if-a-program-exists-from-a-bash-script

# if ! command -v mamba  &> /dev/null
#     then
#         conda install -n base -c conda-forge mamba
# fi

# # Activate base
# eval "$(conda shell.bash hook)"
# conda activate base

# # install and create a snakemake conda environment
# mamba create -c conda-forge -c bioconda -n snakemake snakemake

# # activate the snakemake environment
# eval "$(conda shell.bash hook)"
# conda activate snakemake


#run the pipeline 
 snakemake -p  --jobs 26  --use-conda  --cluster-config \
 /lustre/shared/wfsr-mcfa/projects/internships/luka/viral_metagenomics_pipeline/configs/slurm.yaml \
 --cluster "sbatch  --job-name={cluster.job-name} --cpus-per-task={cluster.cpus-per-task} \
 --qos={cluster.qos} " \
 --rerun-incomplete --keep-incomplete --slurm 


#  --use-conda  --rerun-incomplete  #--cores 16
#--mem-per-cpu {cluster.mem-per-cpu} 
#-n {cluster.ntasks}
#-p {cluster.partition}
#--output {cluster.output} --error {cluster.error} 
#--mem={cluster.mem_mb}
#-n {cluster.ntasks}
# -N {cluster.nodes}
#snakemake --profile
# --mail-type={cluster.mail-type} --mail-user={cluster.mail-user}
# --output={cluster.output} --error={cluster.error}