rule feature_counts:
    input:
        bam = expand("results/mapped/bowtie2/{sample}.bam", sample = SAMPLES),
        gtf = "../resource/ref/genome.gtf",
        fasta = "../resource/ref/genome.fasta"
    output:
        "results/feature_counts/raw_count_matrix_ensemblid.tsv"
    log:
        "logs/feature_counts/feature_counts.log"
    threads:
        config["params"]["general"]["threads"]
    shell:
        """
        featureCounts -T {threads} -t gene -g gene_id -F GTF -a {input.gtf} -O {input.fasta} -J -o {output} {input.bam} 2> {log}
        """