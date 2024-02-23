
process BEDTOOLS_INTERSECT {
    """
    """
    tag { meta.id }
    container 'nciccbr/ccbr_ubuntu_base_20.04:v6'

    input:
    tuple val(meta), path(atac_bed), path(promoter_bed)

    output:
    tuple val(meta), path("intersect.*.bed"),                                emit: intersect
    tuple val(meta), path("count.atac-*.bed"), path("count.promoter-*.bed"), emit: counts

    script:
    def outfile_intersect = "intersect.${atac_bed.baseName}.${promoter_bed.baseName}.bed"
    def outfile_countA = "count.atac-promoter.${atac_bed.baseName}.${promoter_bed.baseName}.bed"
    def outfile_countB = "count.promoter-atac.${promoter_bed.baseName}.${atac_bed.baseName}.bed"
    """
    bedtools intersect \\
        -a ${atac_bed} \\
        -b ${promoter_bed} \\
        -wo \\
        > ${outfile_intersect}
    bedtools intersect \\
        -a ${atac_bed} \\
        -b ${promoter_bed} \\
        -c \\
        > ${outfile_countA}
    bedtools intersect \\
        -a ${promoter_bed} \\
        -b ${atac_bed} \\
        -c \\
        > ${outfile_countB}
    """
}
