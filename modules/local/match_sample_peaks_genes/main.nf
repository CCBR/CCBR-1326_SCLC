process MATCH_SAMPLES_PEAKS_GENES {
    """
    Match atac-seq and rna-seq samples and write to separate files per peak-gene pair
    """
    container "nciccbr/sclc_r-quarto:v0.4.0"

    input:
        tuple path(metadata_infile), path(pdx_infile), path(rna_counts_infile), path(atac_counts_infile), path(peak_gene_infile)

    output:
        path("matched_*.tsv")

    script:
    """
    match_sample_peaks_genes.R ${metadata_infile} ${pdx_infile} ${rna_counts_infile} ${atac_counts_infile} ${peak_gene_infile}
    """

}
