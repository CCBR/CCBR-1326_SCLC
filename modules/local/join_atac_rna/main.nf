
process JOIN_ATAC_RNA {
    """
    Join ATAC promoters with RNAseq DGE
    """
    tag { meta.id }
    container 'nciccbr/consensus_peaks:v1.1'

    input:
    tuple val(meta), path(promoters_logfc), path(rna_logfc)

    output:
    tuple val(meta), path("*.joined_promoters.tsv")

    script:
    outfile = "${rna_logfc.baseName}.joined_promoters.tsv"
    template 'join_atac_rna.R'
}
