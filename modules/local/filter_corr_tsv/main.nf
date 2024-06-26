process FILTER_CORR_TSV {
    """
    Filter gene-peak correlations
    """
    container "nciccbr/sclc_r-quarto:v0.3.3"

    input:
        path(tsv)

    output:
        path("*_filt.tsv")

    script:
    outfile = "${tsv.baseName}_filt.tsv"
    "filter_corr_tsv.R ${tsv} ${outfile}"

    stub:
    """
    touch ${tsv.baseName}_filt.tsv
    """
}
