#runs BLASTn
import os

# """
#     Run BLASTn to search for viral sequences in contigs.

#     Input:
#         contigs (file): Path to the input FASTA file containing contigs.
        
#     Output:
#         blast_search (file): Path to the BLASTn search results in TSV format.

#     Params:
#         db (str): Path to the BLAST database.
#         taxid_db (str): Path to the taxonomic database for BLAST.
#         max_target_seqs (int): Maximum number of target sequences to report.

#     Resources:
#         mem_mb (function): Memory resources allocated for the job in MB.
#         runtime (function): Estimated runtime for the job.

#     Threads:
#         Number of threads to use for BLASTn.

#     Conda:
#         Dependency specification for the Conda environment.

#     Shell:
#         The shell command to execute the BLASTn search.

#     Dependencies:
#         Requires external BLAST  and databases to be configured properly.
    
# """



# rule all:
#     input: 
#         expand(os.path.join(f"{blast_dir}", "bin_size1000", "{bin_id}_blast_out.tsv"), bin_id=BINS)


rule blastn:
    input: 
        contigs = os.path.join(f"{genomad_dir}", "host_removed", "none_checkv_filtered", "{sample}","contigs_summary", "contigs_virus.fna"),
    output: 
        blast_search = os.path.join(f"{blast_dir}", "host_removed", "{sample}", "blastn_search.tsv")

    params:
        db = config["blast_params"]["db"],
        taxid_db = config["blast_params"]["taxid_db"],
        max_target_seqs = config["blast_params"]["max_target_seqs"]


    resources:
        mem_mb = lambda wildcards, attempt: attempt *  config["blast_params"]["mem_mb"],
        runtime =  lambda wildcards, attempt: attempt * config["blast_params"]["runtime"]

    threads:config["blast_params"]["threads"]

    conda: config["conda"]["blast_conda"]

    shell:
        """
       set +u ; export BLASTDB=$BLASTDB:{params.taxid_db}; set -u && blastn \
            -query {input.contigs} \
            -db {params.db} \
            -out {output.blast_search} \
            -max_target_seqs {params.max_target_seqs} \
            -num_threads {threads} \
            -outfmt "6 qseqid staxids sacc sscinames  sblastnames evalue bitscore pident"
        """



