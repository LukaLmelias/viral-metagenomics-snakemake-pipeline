# first concatenate all assembled contigs:
import os



# # include some python scripts
# include: './helpers.py'

# # configs
# configfile: "/lustre/shared/wfsr-mcfa/projects/internships/luka/viral_metagenomics_pipeline/configs/config.yaml"

# metaspades_dir = config['directories']['metaspades_results']
# krakentools_dir = config["directories"]["krakentools_results"]
# vamb_dir = config["directories"]["vamb_results"]

# # extract the path to the raw reads from the config file.
# raw_reads_path  = config["directories"]["raw_reads"]

# SAMPLES,  EXTENSION = glob_wildcards(os.path.join(raw_reads_path, '{sample}_{extension}'))



# rule all:
#     input: 
#         os.path.join(vamb_dir, "allContigsK127.fasta"),
#         os.path.join(vamb_dir, "allContigsK127_sample_ids.tsv"),
#         os.path.join(vamb_dir, "allContigsK127.mmi"),
#         expand(os.path.join(vamb_dir, f"allContigsK127_{{sample}}.bam"), sample=SAMPLES),
#         os.path.join(vamb_dir, "bin_size1000", "clusters.tsv"),

        



# """
#     Concatenate contig files.

#     Parameters:
#     - `input`: Input directory containing  contig files (host-removed in this case but can be any).
#     - `output`: Output concatenated contig file and sample IDs file.
#     - `threads`: Number of threads to use.
#     - `shell`: Python script to concatenate contigs.
#     """


rule concatenate_contigs:    
    input: os.path.join(metaspades_dir, "host_removed")
        

    output: 
        all_contigs_file = os.path.join(vamb_dir,  "allContigsK127.fasta"),   
        sample_ids_file = os.path.join(vamb_dir,  "allContigsK127_sample_ids.tsv")

    threads: 16
    shell:
        """
        python concatenateContigs.py {input} {output.all_contigs_file} {output.sample_ids_file}
        """


# """
#     Create an index for the concatenated contig file using Minimap2.

#     Parameters:
#     - `input`: Input concatenated contig file.
#     - `output`: Output Minimap2 index file.
#     - `params`: Additional parameters (currently empty).
#     - `resources`: Memory and runtime requirements.
#     - `threads`: Number of threads to use.
#     """
rule minimap2Index:    
    input: os.path.join(vamb_dir,  "allContigsK127.fasta"), 

    output: os.path.join(vamb_dir, "allContigsK127.mmi")

    params:

    resources:
        mem_mb = config["minimap2_params"]["mem_mb"],
        runtime = config["minimap2_params"]["runtime"]

    threads: config["minimap2_params"]["threads"]

    conda: config["conda"]["minimap2_conda"]

    shell:
        """
        minimap2 -d {output} {input}
        """

# """
#     Map extracted reads to the concatenated contigs using Minimap2.

#     Parameters:
#     - `input`: Input Minimap2 index file and extracted read files.
#     - `output`: Output SAM file for each sample.
#     - `params`: Additional parameters (currently empty).
#     - `resources`: Memory and runtime requirements.
#     - `threads`: Number of threads to use.
#     - `conda`: Conda environment specification.
#     - `shell`: Command to run Minimap2 for read mapping.
#     """
rule minimap2:    
    input: 
        index = os.path.join(vamb_dir, "allContigsK127.mmi"),
        r1 = os.path.join(f"{krakentools_dir}", 'host_removed_extractedReads', f"{{sample}}_extracted_R1.fasta"), # Expand at the rule all
        r2 = os.path.join(f"{krakentools_dir}", 'host_removed_extractedReads', f"{{sample}}_extracted_R2.fasta")
        

    output: temp(os.path.join(vamb_dir,  f"allContigsK127_{{sample}}.sam"))

    params:

    resources:
        mem_mb = lambda wildcards, attempt: attempt *  config["minimap2_params"]["mem_mb"],
        runtime = lambda wildcards, attempt: attempt * config["minimap2_params"]["runtime"]

    threads: 20

    conda: config["conda"]["minimap2_conda"]

    shell:
        """
        minimap2 -t {threads} -k 15 -w 10 -ax sr {input.index}  {input.r1} {input.r2} -o {output}
        """


# """
#     Convert SAM files to BAM files and sort them.

#     Parameters:
#     - `input`: Input SAM files.
#     - `output`: Output sorted BAM files.
#     - `params`: Additional parameters (currently empty).
#     - `resources`: Memory and runtime requirements.
#     - `threads`: Number of threads to use.
#     - `conda`: Conda environment specification.
#     """
rule samtools:    
    input: os.path.join(vamb_dir,  f"allContigsK127_{{sample}}.sam")


    output: os.path.join(vamb_dir,  f"allContigsK127_{{sample}}.bam")


    params:

    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["samtools_params"]["mem_mb"],
        runtime = lambda wildcards, attempt: attempt * config["samtools_params"]["runtime"]

    threads: config["samtools_params"]["threads"]

    conda: config["conda"]["samtools_conda"]
    shell:
        """
        samtools view  -F 3584 -b --threads {threads} {input} | samtools sort -n --threads {threads} -O BAM -o {output} 
        """

# """
#     Run VAMB to bin contigs based on read coverage.

#     Parameters:
#     - `output`: Output cluster information file.
#     - `params`: Parameters for VAMB, including input BAM files, output directory.
#     - `resources`: Memory and runtime requirements.
#     - `threads`: Number of threads to use.
#     - `conda`: Conda environment specification.
#     """
rule vamb:    
    output: os.path.join(vamb_dir, "bin_size1000", "clusters.tsv")

    params:
        bamfiles = f"{vamb_dir}/*.bam",
        outdir= os.path.join(vamb_dir, "bin_size1000"),
        contigs = os.path.join(vamb_dir, "allContigsK127.fasta"),
        minfasta = config["vamb_params"]["minfasta"],
        min_seq_len = config["vamb_params"]["min_seq_len"],
        initial_batch_size = config["vamb_params"]["initial_batch_size"]
    
    resources:
        mem_mb = config["vamb_params"]["mem_mb"],
        runtime = config["vamb_params"]["runtime"]

    threads: config["vamb_params"]["threads"]
    
    conda: config["conda"]["vamb_conda"]
    
    shell:
        """
      vamb  --outdir {params.outdir} --fasta {params.contigs} -p {threads} \
        --bamfiles {params.bamfiles} \
        -t {params.initial_batch_size} -m {params.min_seq_len} \
        --minfasta {params.minfasta}
        """

# set +u sets off batch strict mode to avoid error on .bai files



# VAMB give this persistent error of file exists even when it doesnt when run with conda:
#File "/lustre/shared/wfsr-mcfa/projects/internships/luka/viral_metagenomics_pipeline/scripts/.snakemake/conda/852adb6864c64b07061d77ce0a48cc76_/lib/python3.7/site-packages/vamb/__main__.py", line 373, in main
    #raise FileExistsError(args.outdir)

# to stop it i hashed out the lines below; unless you have another solution; hash them too before runing vamb
 

 ######################### CHECK INPUT/OUTPUT FILES #####################

    # Outdir does not exist
    # args.outdir = os.path.abspath(args.outdir)
    # if os.path.exists(args.outdir):
    #     raise FileExistsError(args.outdir)








#rule 
#snakemake -n -s vamb.smk   --use-conda --cores 4
#snakemake -n -j 20  --profile slurmProfile -s vamb.smk
