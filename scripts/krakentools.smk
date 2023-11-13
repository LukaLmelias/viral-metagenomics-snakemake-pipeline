
rule krakentools:
    input:
         KRAKEN_FILE = os.path.join(f'{kraken2_dir}', 'results', f'{{sample}}.kraken2_out'),
         SEQ_FILE1 = os.path.join(f'{fastp_dir}', 'qc_reads', f'{{sample}}_R1.qc.fastq.gz'),
         SEQ_FILE2 = os.path.join(f'{fastp_dir}', 'qc_reads', f'{{sample}}_R2.qc.fastq.gz'),
         REPORT_FILE = os.path.join(f'{kraken2_dir}', 'results', f'{{sample}}.kraken2_report'),
         #SEQ_FILE2 = f'{fastp_dir}qc_reads/{{sample}}_R2.qc.fastq.gz',
         
         

    #Documentation:
        #      -k KRAKEN_FILE        Kraken output file to parse
        #   -s SEQ_FILE1, -s1 SEQ_FILE1, -1 SEQ_FILE1, -U SEQ_FILE1
        #                         FASTA/FASTQ File containing the raw sequence letters.
        #   -s2 SEQ_FILE2, -2 SEQ_FILE2
        #                         2nd FASTA/FASTQ File containing the raw sequence letters (paired).
        #   -t TAXID [TAXID ...], --taxid TAXID [TAXID ...]
        #                         Taxonomy ID[s] of reads to extract (space-delimited)


    output:
        OUTPUT_FILE = os.path.join(f"{krakentools_dir}", 'host_removed_extractedReads', f"{{sample}}_extracted_R1.fasta"), # Expand at the rule all
        OUTPUT_FILE2 = os.path.join(f"{krakentools_dir}", 'host_removed_extractedReads', f"{{sample}}_extracted_R2.fasta")


    # Documentation:
        #     -o OUTPUT_FILE, --output OUTPUT_FILE
        #                         Output FASTA/Q file containing the reads and sample IDs
        #      -o2 OUTPUT_FILE2, --output2 OUTPUT_FILE2
        #                         Output FASTA/Q file containig the second pair of reads [required for paired input]    

    params:        
        includeParents  = config['krakentools_params']['includeParents'],
        includeChildren = config['krakentools_params']['includeChildren'],
        TAXID = config["krakentools_params"]["TAXID"]

    #Documentation:
        #     -r REPORT_FILE, --report REPORT_FILE
        #                         Kraken report file. [required only if --include-parents/children is specified]
        #   --include-parents     Include reads classified at parent levels of the specified taxids
        #   --include-children    Include reads classified more specifically than the specified taxids

    conda: config["conda"]["krakentools_conda"]

    resources:
        mincpus = config['krakentools_params']['threads'],
        runtime = config['krakentools_params']['time'],
        mem_mb = config["krakentools_params"]["mem_mb"],

    threads: config['krakentools_params']['threads']

    shell:
        """
        extract_kraken_reads.py -k {input.KRAKEN_FILE} -s {input.SEQ_FILE1} \
        -r {input.REPORT_FILE}  -t {params.TAXID}  --exclude --include-children \
        -o {output.OUTPUT_FILE} -o2 {output.OUTPUT_FILE2} -s2 {input.SEQ_FILE2}

        """
#   snakemake -n  -s krakentools.smk   --use-conda --cores 16

