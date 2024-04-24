process CORRELATE_PEAKS_GENES {
    """
    Correlate ATAC-seq peak counts near genes' TSSs with the genes' RNA-seq counts
    """
    container "nciccbr/sclc_r-quarto:v0.3.3"
    cpus 16
    memory '200 GB'
    //queue 'largemem'
    time '24 h'

    input:
        tuple path(metadata_infile), path(rna_counts_infile), path(atac_counts_infile), path(peak_gene_infile)

    output:
        path("peaks_genes_corr.csv")

    script:
    peak_gene_outfile = "peaks_genes_corr.csv"
    template 'correlate_peaks_genes.R'

    stub:
    """
    touch peaks_genes_corr.csv
    """
}
