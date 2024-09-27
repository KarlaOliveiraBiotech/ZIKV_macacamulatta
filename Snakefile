##### Importing libraries #####
import pandas as pd
from snakemake.utils import min_version
import snakemake_wrapper_utils
 
min_version("7.6.1")
##### Importing libraries #####

##### Loading config file #####
configfile: "../config/config.yaml"
##### Loading config file #####

##### Data structure #####
# resource/data/raw_data/subset/sample_name/sample_file - each sample is inside its own folder

##### Global variable: Getting folders and all prefix for sample filename #####
# Gets samples prefix for Fastp. In case another structure, change accordingly.
FOLDERS, SAMPLES = glob_wildcards("../resource/data/raw_data/subset/{folder}/{sample}_2M.fastq.gz") 



# Gets name of sample files for Bowtie2
#SAMPLES_MAP = glob_wildcards("results/fastp/fq/{sample}.fastp.fq").sample


# Gets name of sample files for featureCounts
#SAMPLES_BAM = glob_wildcards("results/mapped/bowtie/{sample}.bam").sample
# ##### Global variable: Getting folders and all prefix for sample filename #####    

rule all:
    input:
        #expand("results/fastp/fq/{sample}.fastp.fq", sample = SAMPLES)
        #"../resource/ref/genome.fasta"
        #"../resource/ref/transcriptome.fasta"
        #"../resource/ref/genome.gtf"
        #"resource/ref/ncrna.fasta"
        #"results/mapped/bowtie2/{sample}.bam"
        "results/feature_counts/raw_count_matrix_ensemblid.tsv"
        
       

         
# ##### Loading rules #####
include: "rules/fastp.smk"
include: "rules/ref.smk"
include: "rules/bowtie2_build.smk"
include: "rules/bowtie2_align.smk"
include: "rules/feature_counts.smk"



