process HINT_FOOTPRINTING {
    //container "nciccbr/sclc_rgt:v0.1.0"
    module 'rgt'
    cpus 32
    memory '32 GB'
    time '12 h'


    input:
    tuple val(meta), path(bam), path(bai), path(bed)
    path(rgtdata)

    output:
    tuple val(meta), path("${meta.id}_footprints/${meta.id}_footprints/*.bed"),  emit: bed
    tuple val(meta), path("${meta.id}_footprints/${meta.id}_footprints/*.info"), emit: info

    script:
    """
    mkdir ${meta.id}_footprints/
    rgt-hint footprinting --atac-seq --paired-end --organism=hg19 \\
        --output-location=${meta.id}_footprints/ \\
        --output-prefix=${meta.id} \\
        ${bam} ${bed}
    """
}
