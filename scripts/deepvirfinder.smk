import pandas as pd
# # runs DeepVirFinder
# configfile: "/lustre/shared/wfsr-mcfa/projects/internships/luka/viral_metagenomics_pipeline/configs/config.yaml",


# rule all:
#     input:
#         pred = "/lustre/shared/wfsr-mcfa/projects/internships/luka/viral_metagenomics_pipeline/Oesters_results/metaviralspades/allContigs_gt1bp_dvfpred.txt"


rule deepvirfinder:
    input:
        bins = "/lustre/shared/wfsr-mcfa/projects/internships/luka/viral_metagenomics_pipeline/scripts/test/allContigsk127.fasta"

        #params:
            # -i INPUT_FA, --in=INPUT_FA
            #                 input fasta file
        
    output:
        pred = os.path.join(f"{deepvirfinder_dir}", "allContigsk127.fasta._gt1bp_dvfpred.txt" )


        

    params:
        deepVirFinderScript = config['deepvirfinder_params']['deepVirFinderScript'],
        modelsDir = config['deepvirfinder_params']['modelsDir'],
        outDir = config["directories"]["deepvirfinder_results"]


        #params:
            # -c CORE_NUM, --core=CORE_NUM (this is the same as threads)
            #                 number of parallel cores (default 1)
            # -m MODDIR, --mod=MODDIR
            #             model directory (default ./models)

    conda: config["conda"]["deepvirfinder_conda"]

    threads:config["deepvirfinder_params"]["threads"]

    resources:
        mincpus = config['deepvirfinder_params']['threads'],
        runtime = config['deepvirfinder_params']['time'],
        mem_mb = config["deepvirfinder_params"]["mem_mb"],


    shell:
        """
        python {params.deepVirFinderScript} -i {input.bins} -o {params.outDir} -c {threads} -m {params.modelsDir}

        """
# rule filterDeepVirFinder:
#     input:
#         pred = "/lustre/shared/wfsr-mcfa/projects/internships/luka/viral_metagenomics_pipeline/Oesters_results/metaviralspades/before_rr.fasta_gt1bp_dvfpred.txt"
    
#     output:
#         selected = config["directories"]["deepvirfinder_results"] + "read_selected.txt"

#     params:
#         score_cutoff = config["deepvirfinder_params"]["score_cutoff"]
#         pval_cutoff = config["deepvirfinder_params"]["pval_cutoff"]
    
#     run:
#         selected = filterDeepVirFinder(input.pred, params.score_cutoff, params.pval_cutoff)
#         selected.to_csv(config["directories"]["deepvirfinder_results"]+"read_selected.txt")







#commadnline test example
    # python DeepVirFinder/dvf.py  \
    # -i /lustre/shared/wfsr-mcfa/projects/internships/luka/viral_metagenomics_pipeline/scripts/test/metabat2/test_contig.11.fa \
    # -c 1 -o /lustre/shared/wfsr-mcfa/projects/internships/luka/viral_metagenomics_pipeline/scripts/test/deepvirfinder \
    # -m /lustre/shared/wfsr-mcfa/projects/internships/luka/viral_metagenomics_pipeline/tools/DeepVirFinder/models



#   snakemake -n  -s deepvirfinder.smk   --use-conda --cores 16

