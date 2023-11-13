#!/bin/bash


### Tiny script to install a given tool and creates its env

# it uses mamba by default

# it does do not ask for your confirmation after because of -y

echo hi, tool name?:
read tool
mamba create -c conda-forge -c bioconda -n $tool $tool -y

