
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

process BEDTOOLS_INTERSECT_FULL {
    """
    intersect bedfiles and only report intersections
    where 100% of a query in bed_A is overlapped by a record in bed_B
    """
    tag { meta.id }
    container 'nciccbr/ccbr_ubuntu_base_20.04:v6'

    cpus 1
    memory '10 GB'
    time '1h'

    input:
    tuple val(meta), path(bed_A), path(bed_B)

    output:
    tuple val(meta), path("intersect.${bed_A.baseName}.${bed_B}"),   emit: intersect

    script:
    def outfile_intersect = "intersect.${bed_A.baseName}.${bed_B}"
    """
    bedtools intersect \\
        -a ${bed_A} \\
        -b ${bed_B} \\
        -f 1.0 \\
        -wo \\
        | sort | uniq \\
        > ${outfile_intersect}
    """
}

process BEDTOOLS_INTERSECT_WA {
    """
    intersect bedfiles and only report intersections
    where 100% of a query in bed_A is overlapped by a record in bed_B.
    only write original entry in A for each overlap.
    """
    tag { meta.id }
    container 'nciccbr/ccbr_ubuntu_base_20.04:v6'

    input:
    tuple val(meta), path(bed_A), path(bed_B)

    output:
    tuple val(meta), path("intersect.${bed_A.baseName}.${bed_B}"),   emit: intersect

    script:
    def outfile_intersect = "intersect.${bed_A.baseName}.${bed_B}"
    """
    bedtools intersect \\
        -a ${bed_A} \\
        -b ${bed_B} \\
        -f 1.0 \\
        -wa \\
        | sort | uniq \\
        > ${outfile_intersect}
    """
}
