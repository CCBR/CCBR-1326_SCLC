process HINT_FOOTPRINTING {
    //container "nciccbr/sclc_rgt:v0.1.0"

    input:
    tuple val(meta), path(bam), path(bai), path(bed)
    path(rgtdata)

    output:
    tuple val(meta), path("footprints/${meta.id}*")

    script:
    """
    ls $rgtdata
    mkdir footprints
    rgt-hint footprinting --atac-seq --paired-end --organism=hg19 \\
        --output-location=footprints \\
        --output-prefix=${meta.id} \\
        ${bam} ${bed}
    """
}
