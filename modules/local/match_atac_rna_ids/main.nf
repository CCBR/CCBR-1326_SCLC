process MATCH_ATAC_RNA_IDS {
    """
    Match atac-seq and rna-seq samples and write out tidy tsv files
    """
    container "nciccbr/sclc_r-quarto:v0.4.0"
    memory '200 GB'
    time '24 h'

    input:
        tuple path(metadata_infile), path(pdx_infile), path(rna_counts_infile), path(atac_counts_infile)

    output:
        tuple path("matched_sample_ids.tsv"), path("rna_counts_norm_long.tsv"), path("atac_counts_norm_long.tsv")


    script:
    """
    match_atac_rna_sample_ids.R ${metadata_infile} ${pdx_infile} ${rna_counts_infile} ${atac_counts_infile}
    """

}
