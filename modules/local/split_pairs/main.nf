process SPLIT_PAIRS {
    """
    Match motifs with atac & rna counts for each sample
    """
    container "nciccbr/sclc_r-quarto:v0.4.0"
    memory '2 GB'
    time '1 h'
    cpus 8

    input:
        path(infile)

    output:
        path("matches/")

    script:
    """
    split_pairs.R ${infile} ${task.cpus}
    """
}
