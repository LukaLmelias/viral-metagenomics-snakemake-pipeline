#!/bin/bash

# runs my snakemake pipeline with a few configurations for the oesters samples.
snakemake -np -j 36  --profile slurmProfile #--until  tax_metagenome
 #--use-conda --cores 16

#snakemake -p -j 36  --profile slurmProfile  --until tax_metagenome

    # --config "directories"="{results: /lustre/shared/wfsr-mcfa/projects/internships/luka/viral_metagenomics_pipeline/Oesters_results/fastp/, }" 



