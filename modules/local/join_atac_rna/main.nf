
process JOIN_ATAC_RNA {
    """
    Join ATAC promoters with RNAseq DGE
    """
    tag { meta.id }
    container 'nciccbr/consensus_peaks:v1.1'

    input:
    tuple val(meta), path(promoters_logfc), path(rna_logfc)

    output:
    tuple val(meta), path("*.joined_promoters.tsv"), emit: tsv
    tuple val(meta), path("*.log"),                  emit: log

    script:
    outfile = "${rna_logfc.baseName}.joined_promoters.tsv"
    logfile = 'joined_atac_rna.log'
    template 'join_atac_rna.R'
}
