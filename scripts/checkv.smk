import os



downloadDB = config["checkv_params"]["downloadDB"]
if downloadDB:
    rule downloadCheckvDb:
        params:
            checkvDbDir = checkv_dir,        
        threads: 4        
        conda: config["conda"]["checkv_conda"]
        resources:
            mem_mb =  8000,
            runtime = 30
        
        shell:
            """
            checkv download_database {params.checkvDbDir}
            """

# rule all:
#     input:
#         os.path.join(f"{checkv_dir}", "quality_summary.tsv") 
#         #"/lustre/shared/wfsr-mcfa/projects/internships/luka/viral_metagenomics_pipeline/Oesters_results/checkv/quality_summary.tsv"


# """
#     Run CheckV to assess the quality of contigs.

#     Input:
#         contigs (file): Path to the input FASTA file containing contigs.

#     Output:
#         report (file): Path to the CheckV quality report.

#     Params:
#         checkvMode (str): CheckV mode (e.g., 'lineage_wf').
#         checkvDb (str): Path to the CheckV database.
#         outputDir (str): Output directory for CheckV results.

#     Resources:
#         mincpus (int): Minimum number of CPUs required.
#         mem_mb (function): Memory resources allocated for the job in MB.
#         runtime (function): Estimated runtime for the job.

#     Threads:
#         Number of threads to use for CheckV.

#     Conda:
#         Dependency specification for the Conda environment.
#     """"
rule checkv:    
    input: 
        contigs = os.path.join(f"{metaspades_dir}", "host_removed", "{sample}", "contigs.fasta")

    output: 
        report= os.path.join(f"{checkv_dir}", "host_removed", "{sample}", "quality_summary.tsv") 

    params:  
        checkvMode = config["checkv_params"]["checkvMode"],
        checkvDb = config["checkv_params"]["checkvDb"],
        outputDir = os.path.join(f"{checkv_dir}", "host_removed", "{sample}"),


    resources:
        mincpus=config["checkv_params"]["threads"],
        mem_mb= lambda wildcards, attempt: attempt * config["checkv_params"]["mem_mb"],
        runtime=  lambda wildcards, attempt: attempt * config["checkv_params"]["runtime"]

    threads: config["checkv_params"]["threads"]
    conda: config["conda"]["checkv_conda"]

    shell:
        """
        checkv {params.checkvMode} {input.contigs} {params.outputDir} -t {threads} -d {params.checkvDb}
        
        """


# """
#     Select high-quality contig nodes based on CheckV quality assessment.

#     Input:
#         checkv_quality_report (file): Path to the CheckV quality assessment report for a sample.

#     Output:
#         checkv_selected_nodes (file): Path to the file containing selected contig nodes.

#     Params:
#         completeness_cutoff (float): The completeness score cutoff for selecting contig nodes.

#     """
rule checkv_select_nodes:    
    input: 
        checkv_quality_report = os.path.join(f"{checkv_dir}", "host_removed", "{sample}", "quality_summary.tsv") 

    output: 
        checkv_selected_nodes = os.path.join(f"{checkv_dir}", "host_removed", "{sample}", "checkv_medium_to_high_nodes.txt")
    
    params:
        completeness_cutoff = config["checkv_params"]["completeness_cutoff"]    

    
    shell:
        """
        python selected_checkv_nodes.py {input.checkv_quality_report} {params.completeness_cutoff} {output.checkv_selected_nodes} 
        """


# """
#     Filter metagenomic contigs based on CheckV quality assessment.

#     Input:
#         checkv_selected_nodes (file): Path to the selected contig nodes file from CheckV.
#         contigs (file): Path to the input FASTA file containing metagenomic contigs.

#     Output:
#         selected_contigs (file): Path to the filtered contigs based on CheckV quality.

#     Conda:
#         Dependency specification for the Conda environment.
#     """
rule filterContigs_on_checkv_quality:    
    input:
        checkv_selected_nodes = os.path.join(f"{checkv_dir}", "host_removed", "{sample}", "checkv_medium_to_high_nodes.txt"),
        contigs = os.path.join(f"{metaspades_dir}","host_removed", "{sample}", "contigs.fasta"),
    output:
        selected_contigs = os.path.join(f"{checkv_dir}", "host_removed", "{sample}", "checkv_medium_to_high_contigs.fasta"),
    
    conda: config["conda"]["seqtk_conda"]

    shell:
        """
        seqtk subseq {input.contigs} {input.checkv_selected_nodes} > {output.selected_contigs}
        """



#snakemake -n -j 20  --profile slurmProfile  -s checkv.smk   





