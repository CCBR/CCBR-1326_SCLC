
process BEDTOOLS_INTERSECT {
    """
    """
    tag { meta.id }
    container 'nciccbr/ccbr_ubuntu_base_20.04:v6'

    input:
    tuple val(meta), path(bedA), path(bedB)

    output:
    tuple val(meta), path("intersect.*.bed")

    script:
    def outfile = "intersect.${bedA.baseName}.${bedB.baseName}.bed"
    """
    bedtools intersect \\
        -a ${bedA} \\
        -b ${bedB} \\
        -wo \\
        > ${outfile}
    """
}
