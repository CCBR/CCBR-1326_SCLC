
process JOIN_PEAKS_PROMOTERS {
    """
    Join atac peaks overlapping promoters with logFC data
    """
    tag { meta.id }
    container 'nciccbr/consensus_peaks:v1.1'

    input:
    tuple val(meta), path(atac_promoters), path(atac_logfc), path(rna_dge)

    output:
    tuple val(meta), path("*.tsv")

    script:
    outfile = "${atac_promoters.baseName}.tsv"
    template 'join_peaks_promoters.R'
}