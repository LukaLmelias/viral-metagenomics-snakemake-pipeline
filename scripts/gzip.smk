
# rule to gzip files to reduce size

rule gzip:
    input: 
        path = '/lustre/shared/wfsr-mcfa/projects/internships/luka/viral_metagenomics_pipeline/results/fastp/qc_reads' # would be nice to make it in a config

    output:
        r1 = f'{fastp_dir}qc_reads/{{sample}}_R1.qc.fastq.gz',
        r2 = f'{fastp_dir}qc_reads/{{sample}}_R2.qc.fastq.gz'
    shell:
        """
        gzip -rv {input.path}
        """