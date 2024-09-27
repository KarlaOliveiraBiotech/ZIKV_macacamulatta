rule gene2symbol:
    input:
        counts = "../results/feature_counts/raw_count_matrix_ensemblid.tsv"
    output:
        raw = 