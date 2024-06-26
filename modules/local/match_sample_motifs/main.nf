process MATCH_SAMPLE_MOTIFS {
    """
    Match motifs with atac & rna counts for each sample
    """
    container "nciccbr/sclc_r-quarto:v0.4.0"
    memory '200 GB'
    time '12 h'
    cpus 8

    input:
        tuple path(motifs_peaks), path(matched_ids), path(rna_counts), path(atac_counts)

    output:
        path("matched_${motifs_peaks.baseName}.tsv")

    script:
    """
    match_sample_motifs.R ${matched_ids} ${rna_counts} ${atac_counts} ${motifs_peaks} matched_${motifs_peaks.baseName}.tsv ${task.cpus}
    """

}
