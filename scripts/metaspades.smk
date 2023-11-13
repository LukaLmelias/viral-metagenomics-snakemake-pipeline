# """
#     Run MetaSPAdes for metagenomic assembly.

#     Parameters:
#     - `r1` and `r2`: Input paired-end FASTA files.
#     - `out_dir`: Output directory for assembled contigs.

#     Conda environment: specified in the config file.

#     Resources:
#     - `mincpus`: Minimum CPU cores required.
#     - `runtime`: Expected runtime for the rule.
#     - `mem_mb`: Memory in megabytes required.

#     Output:
#     - `out_dir`: Output directory containing assembled contigs.
#     """

rule metaspades:    
    input:
        r1 = os.path.join(f"{krakentools_dir}", 'host_removed_extractedReads', f"{{sample}}_extracted_R1.fasta"),
        r2 = os.path.join(f"{krakentools_dir}", 'host_removed_extractedReads', f"{{sample}}_extracted_R2.fasta")


    output:
        out_dir =  os.path.join(f"{metaspades_dir}","host_removed", "{sample}", "contigs.fasta")
        
        
        

    params:
        out_dir =  os.path.join(f"{metaspades_dir}", "host_removed", "{sample}"),
        checkpoint = 'all',
        kmer = config["metaviraspades_params"]["kmer"]
         # ["all" | "last"]
        #out_dir = config['directories']['metaviralspades_results'],
        # trusted_contigs = '{sample}-trusted-contigs' ,# not sure if it should be absolute path;
        # untrusted_contigs = '{sample}-untrusted-contigs',
        # assembly_graph = '{sample}-assembly-graph'

    conda: config['conda']['metaviralspades_conda']

    threads: config['metaviraspades_params']['threads']

    resources:
        mincpus = config['metaviraspades_params']['threads'],
        runtime = lambda wildcards, attempt: attempt * config['metaviraspades_params']['time'],
        #mem_per_cpu = config["metaviraspades_params"]["mem_per_cpu"],
        mem_mb= lambda wildcards, attempt: attempt * config["metaviraspades_params"]["mem_mb"]
        
    shell:
        """
        metaspades.py  -t {threads} -1 {input.r1} -2 {input.r2} \
         --checkpoints {params.checkpoint} -k {params.kmer} --only-assembler -o {params.out_dir} 
        """

    