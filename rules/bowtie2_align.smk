##### Align fastq samples against reference genome #####

rule bowtie2_align:
    input:
        sample = ["results/fastp/fq/{sample}.fastp.fq"],
        idx = multiext(
            "../resource/ref/index/mmul_10",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    params:
        basename = "../resource/ref/index/mmul_10"
    output: 
        "results/mapped/bowtie2/{sample}.bam"
    log: 
        "logs/bowtie2/{sample}.log"
    threads: 
        config["params"]["general"]["threads"]  # Use at least two threads
    # wrapper: config["wrapper_version"] + "/bio/bowtie2/align"
    shell:
        """
        bowtie2 -p {threads} -x {params.basename} -U {input.sample} 2> {log} | samtools view -bS - > {output}
        """