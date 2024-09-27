rule bowtie2_build:
    input: 
        "../resource/ref/genome.fasta",
    params:
        basename = "../resource/ref/index/mmul_10"
    output:
        multiext(
            "../resource/ref/index/mmul_10",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
        # This structure is pivotal. Without it, rule will not run.
    log: 
        "logs/bowtie2_build/build.log",
    threads: 
        config["params"]["general"]["threads"]
    shell:
        """
        bowtie2-build --threads {threads} {input} {params.basename} 2> {log}
        """
    #wrapper: config["wrapper_version"] + "/bio/bowtie2/build"