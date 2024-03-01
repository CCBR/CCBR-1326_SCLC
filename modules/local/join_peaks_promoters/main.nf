
process JOIN_PEAKS_PROMOTERS {
    """
    Join atac peaks overlapping promoters with logFC data
    """
    tag { meta.id }
    container 'nciccbr/consensus_peaks:v1.1'

    input:
    tuple val(meta), path(atac_promoters), path(atac_logfc), path(rna_logfc)

    output:
    tuple val(meta), path("*.atac-promoters.tsv"), path(rna_logfc)

    script:
    outfile = "${atac_promoters.baseName}.atac-promoters.tsv"
    template 'join_peaks_promoters.R'
}
