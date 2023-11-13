

rule haploflow:
    input:
        "/lustre/shared/wfsr-mcfa/projects/internships/luka/viral_metagenomics_pipeline/Oesters_results/fastp/qc_reads/combined__01_ACACGATC-ATGGTATT_R1.qc.fastq.gz"
        
    output:
        out_dir = directory('/lustre/shared/wfsr-mcfa/projects/internships/luka/viral_metagenomics_pipeline/Oesters_results/haploflow/')
    resources:
        mincpus=config['haploflow_params']['threads'],
        runtime=config['haploflow_params']['time'],
        mem_per_cpu=config["haploflow_params"]["mem-per-cpu"],
        #mem_mb=config["haploflow_params"]["mem_mb"],
        cpus_per_task=config["haploflow_params"]["cpus_per_task"]
    log:
        '/lustre/shared/wfsr-mcfa/projects/internships/luka/viral_metagenomics_pipeline/Oesters_results/haploflow/log'
    conda: config["conda"]["haploflow_conda"]
    shell:
        """
        /lustre/shared/wfsr-mcfa/projects/internships/luka/viral_metagenomics_pipeline/tools/Haploflow/haploflow \
        --read-file {input} --out {output.out_dir} --log {log}
        """