##### Gets all reference information from ENSEMBLE for the specie written in configfile #####


rule genome_ref:
    output:
       "../resource/ref/genome.fasta",
    log:
        "logs/ref/genome_ref.log",
    params:
        species = config["ref"]["species"],
        build = config["ref"]["build"],
        datatype = "dna",
        id_type = "toplevel",
        release = config["ref"]["release"],
    message: "Downloading the genome {params.species} ({params.build}) and release {params.release}"
    wrapper:
        config["wrapper_version"] + "/bio/reference/ensembl-sequence"


rule transcriptome_ref:
    output:
        "../resource/ref/transcriptome.fasta"
    log:
        "logs/ref/transcriptome_ref.log"
    params:
        species = config["ref"]["species"],
        datatype = "cdna",
        build = config["ref"]["build"],
        release = config["ref"]["release"]
    message:
        "Downloading the transcriptome {params.species} ({params.build}) and release {params.release}"
    wrapper: config["wrapper_version"] + "/bio/reference/ensembl-sequence"


rule annotation_ref:
    output:
        "../resource/ref/genome.gtf"
    log:
        "logs/ref/annotation_ref.log"
    params:
        species = config["ref"]["species"],
        datatype = "gtf",
        build = config["ref"]["build"],
        release = config["ref"]["release"]
    message:
        "Downloading the annotation {params.species} ({params.build}) and release {params.release}"
    wrapper: config["wrapper_version"] + "/bio/reference/ensembl-annotation"  


rule ncrna_ref:
    output:
        "../resource/ref/ncrna.fasta"
    log:
        "logs/ref/ncrna_ref.log"
    params:
        species = config["ref"]["species"],
        datatype = "dna",
        build = config["ref"]["build"],
        release = config["ref"]["release"]
    message:
        "Downloading the ncRNA {params.species} ({params.build}) and release {params.release}"
    wrapper: config["wrapper_version"] + "/bio/reference/ensembl-sequence"