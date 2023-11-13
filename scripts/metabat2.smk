
# """
#     Run Metabat2 binning tool on contigs.

#     Parameters:
#     - `contig`: Input contig file.
#     - `outFile`: Output directory for generated bins.
#     - `minContig`: Minimum size of a contig for binning (default: 2500).
#     - `minBinSize`: Minimum size of a bin as the output (default: 200000).
#     - `percentGoodContig`: Percentage of 'good' contigs considered for binning (default: 95).
#     - `threads`: Number of CPU threads to use (default: 8).

#     Conda environment: specified in the config file.

#     Resources:
#     - `mincpus`: Minimum CPU cores required.
#     - `runtime`: Expected runtime for the rule.
#     - `mem_mb`: Memory in megabytes required.

#     Output:
#     - `outFile`: Output directory containing bins and unbinned contigs.
#     """    

rule metabat2:    
    input:
        contig =os.path.join(f"{genomad_dir}", "host_removed", "none_checkv_filtered", "{sample}","contigs_summary", "contigs_virus.fna"),
        
    output:
        outFile = directory(os.path.join(f"{metabat2_dir}", "host_removed_checkv_filtered","{sample}"))
        
        
        #params documentation:        
            # -o [ --outFile ] arg              Base file name and path for each bin
            # -a [ --abdFile ] arg              A file having mean and variance of base coverage depth
            # --unbinned                        Generate [outFile].unbinned.fa file for unbinned contigs



    params:
        minContig = config['metabat2_params']['MinContig'],
        minBinSize = config['metabat2_params']['minBinSize'], # both minBin size and percentGoodContig affetct the resulting bins; 
        percentGoodContig = config['metabat2_params']['percentGoodContig'],
        outFile = os.path.join(f"{metabat2_dir}", "host_removed_checkv_filtered","{sample}","{sample}") 

        #params documentations:
            #-m [ --minContig ] arg (=2500)    Minimum size of a contig for binning (should be >=1500).
            #-s [ --minClsSize ] arg (=200000) Minimum size of a bin as the output.
            # --maxP arg (=95)                  Percentage of 'good' contigs considered for binning decided by connection
            #                         among contigs. The greater, the more sensitive.


    conda: config['conda']['metabat2_conda']


    resources:
        mincpus = config['metabat2_params']['threads'],
        runtime = lambda wildcards, attempt: attempt *config['metabat2_params']['time'],
        #mem_mb_per_cpu = config["metabat2_params"]["mem_per_cpu"],
        mem_mb= lambda wildcards, attempt: attempt * config["metabat2_params"]["mem_mb"]
        

    threads: config['metabat2_params']['threads'],


    shell:
        """        
        metabat2  -i {input.contig} -o {params.outFile} -t {threads}  -s {params.minBinSize} \
        -m {params.minContig}  --maxP {params.percentGoodContig} --unbinned 
        """

#   snakemake -n  -s metabat2.smk   --use-conda --cores 8