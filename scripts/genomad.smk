
# """
#         Download the Genomad database if specified in the configuration.

#         Input:
#             genomadDbDir (directory): Output directory for Genomad database download.

#         Output:
#             db (directory): Directory containing the downloaded Genomad database.
#         """

downloadDB = config["genomad_params"]["downloadDB"]

if downloadDB:
    rule downloadgenomadDb:        
        input:
            genomadDbDir = genomad_dir
        output:
            db = directory(os.path.join(genomad_dir, 'genomad_db'))
        
        threads: 8
        
        conda: config["conda"]["genomad_conda"]

        resources:
            mem_mb =  8000, 
            runtime = 30
        
        shell:
            """
            genomad download-database {input.genomadDbDir}
            """

# """
#     Run Genomad to predict viruses in contigs.

#     Input:
#         contigs (file): Path to the input contigs file.
#         genomadDb (directory): Path to the Genomad database directory.

#     Output:
#         viruses (file): Path to the output file containing predicted viruses.

#     Params:
#         genomadMode (str): Genomad mode specified in the configuration.
#         outputDir (directory): Output directory for Genomad results.
#     """
rule genomad:    
    input: 
        contigs = os.path.join(f"{metaspades_dir}","host_removed", "{sample}", "contigs.fasta"),
        genomadDb = os.path.join(genomad_dir, 'genomad_db')

    output: 
        viruses = os.path.join(f"{genomad_dir}", "host_removed", "none_checkv_filtered", "{sample}","contigs_summary", "contigs_virus.fna"),

    params: 
        genomadMode = config["genomad_params"]["genomadMode"],        
        outputDir = os.path.join(f"{genomad_dir}", "host_removed", "none_checkv_filtered", "{sample}")

    resources:
        mincpus=config["genomad_params"]["threads"],
        mem_mb= lambda wildcards, attempt: attempt * config["genomad_params"]["mem_mb"],
        runtime= lambda wildcards, attempt: attempt * config["genomad_params"]["runtime"]

    threads: config["genomad_params"]["threads"]

    conda: config["conda"]["genomad_conda"]


    shell:
        """
        genomad {params.genomadMode} -t {threads} {input.contigs} {params.outputDir} {input.genomadDb}
        """
#snakemake -p -j 1  -s genomad.smk  --use-conda --cores 16
