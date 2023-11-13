# runs kraken2
import os


# """
#     Run Kraken2 to classify sequences in paired-end FASTQ files.

#     Input:
#         r1 (file): Path to the input R1 FASTQ file.
#         r2 (file): Path to the input R2 FASTQ file.

#     Output:
#         report (file): Path to the Kraken2 report file.
#         Krake_file (file): Path to the Kraken2 output file.

#     Params:
#         db (str): Path to the Kraken2 database.
#         threads (int): Number of threads to use for Kraken2 execution.
#     """
rule kraken2:
    input:
        r1 = os.path.join(f'{fastp_dir}', 'qc_reads', f'{{sample}}_R1.qc.fastq.gz'),
        r2 = os.path.join(f'{fastp_dir}', 'qc_reads', f'{{sample}}_R2.qc.fastq.gz'),

    output: 
        report = os.path.join(f'{kraken2_dir}', 'results', f'{{sample}}.kraken2_report'),
        Krake_file = os.path.join(f'{kraken2_dir}', 'results', f'{{sample}}.kraken2_out'),
             
        

        


        #Documentation:
            # --unclassified-out FILENAME
            #                   Print unclassified sequences to filename
            # --classified-out FILENAME
            #                   Print classified sequences to filename



    params:
        db = config["kraken2_params"]['db'],
        threads = config["kraken2_params"]['threads'],
        
    threads: config["kraken2_params"]['threads']

    log: os.path.join(f'{kraken2_dir}', 'logs', f'{{sample}}.log')

    benchmark: os.path.join(f'{kraken2_dir}', 'benchmark', f'{{sample}}.tsv')

    resources:
        mincpus = config['kraken2_params']['threads'],
        runtime= lambda wildcards, attempt: attempt * config['kraken2_params']['time'],
        #mem_per_cpu = config["kraken2_params"]["mem_per_cpu"],
        mem_mb=lambda wildcards, attempt: attempt * config['kraken2_params']['mem_mb'],
          
    conda: config["conda"]["kraken2_conda"]

    shell:
        """
        time kraken2 --db {params.db} \
        --threads {threads} \
        --use-names   \
        --report {output.report} --output {output.Krake_file}\
        --paired {input.r1} {input.r2} \
        --gzip-compressed  &>> {log}
        """
#&>> {log}


rule kraken2_multiqc:       
    output:
        multiqc_data = directory(os.path.join(f'{kraken2_dir}' , 'multiqc_reports', 'multiqc_data'))
    conda: config['conda']['multiqc_conda']
    
    params:
        out_folder = os.path.join(f'{kraken2_dir}' , 'multiqc_reports', 'multiqc_data'),
        input_dir = os.path.join(f'{kraken2_dir}', 'results')
    message:
        "Summary of Kraken reports"
    shell:
        "multiqc {params.input_dir}  -o {params.out_folder} "

