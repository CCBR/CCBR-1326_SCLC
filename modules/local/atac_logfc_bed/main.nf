
process ATAC_LOGFC_BED {
    """
    Convert ATAC LogFC file to BED6 file
    """
    tag { meta.id }
    container 'nciccbr/consensus_peaks:v1.1'

    input:
    tuple val(meta), path(atac_logfc)

    output:
    tuple val(meta), path("*.bed")

    script:
    infile = atac_logfc
    outfile = "${atac_logfc.baseName}.bed"
    template 'atac_logfc_bed.R'
}