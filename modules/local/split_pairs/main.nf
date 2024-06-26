process SPLIT_PAIRS {
    """
    Match motifs with atac & rna counts for each sample
    """
    container "nciccbr/sclc_r-quarto:v0.4.0"
    memory '20 GB'
    time '12 h'
    cpus 12

    input:
        path(infile)

    output:
        path("matches/")

    script:
    """
    split_pairs.R ${infile} ${task.cpus}
    """

}
