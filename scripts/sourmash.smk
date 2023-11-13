

# #######################################################################################################################
# ######## SOURMASH Pipeline (adopted from: https://github.com/PacificBiosciences/pb-metagenomics-tools)
# ######################################################################################################################


onstart:
    print("------------------------------")
    print("sourmash taxonomic profiling workflow")
    print("------------------------------")


onsuccess:
    print("\n--- Workflow executed successfully! ---\n")

onerror:
    print("Alas!\n")


# rule all:
#     input:
#         expand(os.path.join(out_dir, '2-gather', '{sample}.k{ks}.gather.with-lineages.csv'),sample=SAMPLES, ks=ksize),
#         expand(os.path.join(out_dir, '3-taxprofile', '{sample}.k{ks}.gather.genbank.kreport.txt'), sample=SAMPLES, ks=ksize),


# """
#     Sketch sequences from the NCBI database for Sourmash analysis.

#     Parameters:
#     - `input`: Input NCBI database FASTA file.
#     - `output`: Output Sourmash signature file.
#     - `threads`: Number of threads to use.
#     - `resources`: Memory and runtime requirements.
#     - `log`: Log file for recording the execution.
#     - `benchmark`: Benchmark file for performance tracking.
#     - `conda`: Conda environment specification.
#     """    
rule sourmash_sketch_ncbi:    
    input: os.path.join(sourmash_params['ncbi_db'], 'AllNucleotide.fa') #downloaded from: https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/AllNucleotide/

    output: os.path.join(sourmash_params['ncbi_db'], 'ncbi-all-viruses-2023-08-k{ks}.sig')

    threads: sourmash_params['threads']
    resources:
        mem_mb=lambda wildcards, attempt: attempt *100000,
        runtime=lambda wildcards, attempt: attempt *360,
    log: os.path.join(logs_dir, "sketch", "sketch-ncbi-all-viruses-2023-08-k{ks}.log")

    benchmark: os.path.join(benchmarks_dir, "sketch", "sketch-ncbi-2023-08-k{ks}.benchmark")

    conda: sourmash_params['sourmash_conda']# "envs/sourmash.yml"    


    shell:
        """
        sourmash sketch dna {input} -p k={wildcards.ks},dna,scaled=1000,abund --singleton \
                                    -o {output} 2> {log}
        """
# """
#     Index the NCBI Sourmash signature file.

#     Parameters:
#     - `input`: Input Sourmash signature file.
#     - `output`: Output Sourmash index file.
#     - `threads`: Number of threads to use.
#     - `resources`: Memory and runtime requirements.
#     - `log`: Log file for recording the execution.
#     - `conda`: Conda environment specification.
#     - `params`: Output base for indexing.
#     """
rule sourmash_index_ncbi: 
    input: os.path.join(sourmash_params['ncbi_db'], 'ncbi-all-viruses-2023-08-k{ks}.sig')
    output: os.path.join(sourmash_params['ncbi_db'], 'ncbi-all-viruses-2023-08-k{ks}.sbt.zip')
    threads: sourmash_params['threads']
    resources:
        mem_mb=lambda wildcards, attempt: attempt *300000,
        runtime=10200,
    log: os.path.join(logs_dir, "sketch", "index-ncbi-all-viruses-2023-08-k{ks}.log")

    conda: sourmash_params['sourmash_conda']# "envs/sourmash.yml"    

    params:
        outputbase = "ncbi-all-viruses-2023-08-k{ks}"

    shell:
        """
        sourmash index --dna -k {wildcards.ks} --scaled 1000 {params.outputbase} {input} 2> {log}
        """



# """
#     Sketch sequences from viral contigs for Sourmash analysis.

#     Parameters:
#     - `input`: Input viral contigs FASTA file.
#     - `output`: Output Sourmash signature file.
#     - `threads`: Number of threads to use.
#     - `resources`: Memory and runtime requirements.
#     - `log`: Log file for recording the execution.
#     - `benchmark`: Benchmark file for performance tracking.
#     - `conda`: Conda environment specification.
#     """
rule sourmash_sketch_dna:    
    input: os.path.join(f"{genomad_dir}", "host_removed",  "{sample}","checkv_medium_to_high_contigs_summary", "checkv_medium_to_high_contigs_virus.fna"),
    output:
        os.path.join(out_dir, "1-sketch", "{sample}.dna.sig.zip")
    threads: sourmash_params['threads']
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=lambda wildcards, attempt: attempt *60,
    log: os.path.join(logs_dir, "sketch", "{sample}.sketch_dna.log")
    benchmark: os.path.join(benchmarks_dir, "sketch", "{sample}.sketch_dna.benchmark")
    conda: sourmash_params['sourmash_conda']# "envs/sourmash.yml"
    shell:
        """
        sourmash sketch dna {input} -p k=21,k=31,k=51,dna,scaled=100,abund \
                                    --name {wildcards.sample} -o {output} 2> {log}
        """

# """
#     Perform Sourmash database search and gather results.

#     Parameters:
#     - `input`: Input Sourmash query signature and databases.
#     - `output`: Output gather results in CSV and text formats.
#     - `params`: Threshold for base pairs (default '0').
#     - `threads`: Number of threads to use.
#     - `resources`: Memory and runtime requirements.
#     - `log`: Log file for recording the execution.
#     - `benchmark`: Benchmark file for performance tracking.
#     - `conda`: Conda environment specification.
#     """
rule sourmash_gather:    
    input:
        query=os.path.join(out_dir, "1-sketch", "{sample}.dna.sig.zip"),
        databases = lambda w: search_databases[f"k{w.ksize}"],
    output:
        gather_csv=os.path.join(out_dir, '2-gather', '{sample}.k{ksize}.gather.csv'),
        gather_txt=os.path.join(out_dir, '2-gather', '{sample}.k{ksize}.gather.txt'),
    params:
        threshold_bp = sourmash_params.get('threshold_bp', '0'),
    threads: sourmash_params['threads']
    resources:
        mem_mb=lambda wildcards, attempt: attempt *500000,
        runtime=lambda wildcards, attempt: attempt *360,
    log: os.path.join(logs_dir, "gather", "{sample}.k{ksize}.gather.log")
    benchmark: os.path.join(benchmarks_dir, "gather", "{sample}.k{ksize}.gather.benchmark")
    conda: sourmash_params['sourmash_conda'] #"envs/sourmash.yml"
    shell:
        # touch output to let workflow continue in cases where 0 results are found
        """
        echo "DB(s): {input.databases}"
        echo "DB(s): {input.databases}" > {log}

        sourmash gather {input.query} {input.databases} --dna --ksize {wildcards.ksize} \
                 --threshold-bp {params.threshold_bp} \
                 -o {output.gather_csv} > {output.gather_txt} 2>> {log}
        
        touch {output.gather_txt}
        touch {output.gather_csv}
        """
# """
#     Perform taxonomic annotation using Sourmash gather results.

#     Parameters:
#     - `input`: Input gather results and lineage database files.
#     - `output`: Output taxonomic annotation results in various formats (krona, CSV, kreport).
#     - `threads`: Number of threads to use.
#     - `resources`: Memory and runtime requirements.
#     - `params`: Output directory and base name for the annotation files.
#     - `conda`: Conda environment specification.
#     - `log`: Log file for recording the execution.
#     - `benchmark`: Benchmark file for performance tracking.
# """
rule tax_metagenome:    
    input:
        gather = os.path.join(out_dir, '2-gather', '{sample}.k{ksize}.gather.csv'),
        lineages = sourmash_params['database_lineage_files'],
    output:
        os.path.join(out_dir, '3-taxprofile', '{sample}.k{ksize}.gather.genbank.krona.tsv'),
        os.path.join(out_dir, '3-taxprofile', '{sample}.k{ksize}.gather.genbank.summarized.csv'),
        os.path.join(out_dir, '3-taxprofile', '{sample}.k{ksize}.gather.genbank.kreport.txt'),
    threads: sourmash_params['threads']
    resources:
        mem_mb=lambda wildcards, attempt: attempt *100000,
        runtime=lambda wildcards, attempt: attempt *360,
    params:
        outd= lambda w: os.path.join(out_dir, f'3-taxprofile'),
        out_base= lambda w: f'{w.sample}.k{w.ksize}.gather.genbank',
    conda: sourmash_params['sourmash_conda'] # "envs/sourmash.yml"
    log: os.path.join(logs_dir, "tax_metagenome", "{sample}.k{ksize}.tax_metagenome.log")
    benchmark: os.path.join(benchmarks_dir, "tax_metagenome", "{sample}.k{ksize}.tax_metagenome.benchmark")
    shell:
        """
        mkdir -p {params.outd}
        sourmash tax metagenome -g {input.gather} -t {input.lineages} -o {params.out_base} \
                                --output-dir {params.outd} --output-format krona csv_summary kreport \
                                --rank species 2> {log}
        """


# """
# Annotate gather results with taxonomic information.

# Parameters:
# - `input`: Input gather results and lineage database files.
# - `output`: Output annotated gather results with taxonomic information.
# - `threads`: Number of threads to use.
# - `resources`: Memory and runtime requirements.
# - `params`: Output directory for annotated results.
# - `conda`: Conda environment specification.
# - `log`: Log file for recording the execution.
# - `benchmark`: Benchmark file for performance tracking.
# """

rule tax_annotate:    
    input:
        gather = os.path.join(out_dir, '2-gather', '{sample}.k{ksize}.gather.csv'),
        lineages = sourmash_params['database_lineage_files'],
    output:
        os.path.join(out_dir, '2-gather', '{sample}.k{ksize}.gather.with-lineages.csv'),
    threads: sourmash_params['threads']
    resources:
        mem_mb=lambda wildcards, attempt: attempt *100000,
        time=lambda wildcards, attempt: attempt *240,
    params:
        outd= lambda w: os.path.join(out_dir, f'2-gather'),
    conda: sourmash_params['sourmash_conda'] #  "envs/sourmash.yml"
    log: os.path.join(logs_dir, "tax_annotate", "{sample}.k{ksize}.tax_annotate.log")
    benchmark: os.path.join(benchmarks_dir, "tax_annotate", "{sample}.k{ksize}.tax_annotate.benchmark")
    shell:
        """
        mkdir -p {params.outd}
        sourmash tax annotate -g {input.gather} -t {input.lineages} -o {params.outd} 2> {log}
        """

#snakemake -n -j 1  -s sourmash.smk  --use-conda --cores 16