# rule for checking the quality of the assembly:

rule metaquast:
    input:
        contig = os.path.join(f"{metaspades_dir}", "{sample}", "contigs.fasta")  
        

    output:
        report = os.path.join(f"{metaquast_dir}", "{sample}","report.html") # THIS IS WHAT GOES TO RULE ALL


    params:
        out_dir = os.path.join(f"{metaquast_dir}", "{sample}"),
        db = config['metaquast_params']['reference_db'],
        maxRefNum = config['metaquast_params']['maxRefNum'],
        minContig = config['metaquast_params']['minContig'],
        #assembly = config['metaquast_params']['assembly'] not necessary

    threads: config['metaquast_params']['threads']

    resources:
        mincpus=config['metaquast_params']['threads'],
        time=config['metaquast_params']['time'],
        mem_mb = config['metaquast_params']['mem_mb'],

    conda: config['conda']['metaquast_conda']

    shell:
        """
        metaquast.py --threads {threads} --output-dir {params.out_dir} \
        --max-ref-number  {params.maxRefNum} --min-contig {params.minContig}\
        --blast-db {params.db} {input.contig}
        """

        #{params.assembly}

     #   snakemake -n  -s metaquast.smk   --use-conda --cores 8

 # --max-ref-number      