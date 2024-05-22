process CORRELATE_PAIR {
    tag { gene }
    """
    Correlate ATAC-seq peak counts near genes' TSSs with the genes' RNA-seq counts
    """
    container "nciccbr/sclc_r-quarto:v0.3.3"
    cpus 12
    memory '20 G'

    input:
        tuple val(gene), path(tsvs)

    output:
        path("*_corr.tsv")

    script:
    "correlation.R ${gene} ${tsvs.join(',')} ${task.cpus}"

    stub:
    """
    touch ${tsvs[0].baseName}_corr.tsv
    """
}
