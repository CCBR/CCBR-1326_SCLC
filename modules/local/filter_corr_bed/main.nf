process FILTER_CORR_BED {
    """
    Filter gene-peak correlations and output in BED format
    """
    container "nciccbr/sclc_r-quarto:v0.3.3"

    input:
        path(tsv)

    output:
        path("*_filt.bed")

    script:
    outfile = "${tsv.baseName}_filt.bed"
    "filter_corr_bed.R ${tsv} ${outfile}"

    stub:
    """
    touch ${tsv.baseName}_filt.bed
    """
}
