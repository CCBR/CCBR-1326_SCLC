process CORRELATE_PAIR {
    """
    Correlate ATAC-seq peak counts near genes' TSSs with the genes' RNA-seq counts
    """
    container "nciccbr/sclc_r-quarto:v0.3.3"

    input:
        path(tsv)

    output:
        path("*_corr.tsv")

    script:
    peak_gene_outfile = "${tsv.baseName}_corr.tsv"
    template 'correlation.R'

    stub:
    """
    touch ${tsv.baseName}_corr.tsv
    """
}
