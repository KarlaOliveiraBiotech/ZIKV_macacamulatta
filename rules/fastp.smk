##### Runs fastp using fastq samples ##### 
# Stores fastq, html and json reports #


rule fastp:
    input:
        "../resource/data/raw_data/subset/{sample}/{sample}_2M.fastq.gz"
    output:
        fq = "results/fastp/fq/{sample}.fastp.fq",
        html = "results/fastp/html/{sample}.html",
        json = "results/fastp/json/{sample}.json"
    log:
        log = "logs/fastp/{sample}_fastp.log"
    shell:
        "fastp -i {input} -o {output.fq} -h {output.html} -j {output.json} 2> {log}"