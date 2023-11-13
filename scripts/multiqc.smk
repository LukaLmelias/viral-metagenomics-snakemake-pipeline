
rule multiQC:
    input:
        input_dir = f'{fastp_dir}' + 'jsons/'
    output:
        fastp_html = f'{fastp_dir}' + 'multiqc_reports/multiqc_report.html'
    conda: config['conda']['multiqc_conda']
    log:
        f'{fastp_dir}' + 'multiqc_reports/fastp_multiqc.log'
    params:
        out_folder = f'{fastp_dir}' + 'multiqc_reports/'
    message:
        "Summary of fastp report across samples."
    shell:
        "multiqc {input.input_dir}  -o {params.out_folder} &>> {log}"

