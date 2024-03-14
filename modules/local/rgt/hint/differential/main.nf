process MOTIFANALYSIS_MATCHING {
    //container "nciccbr/sclc_rgt:v0.1.0"
    module 'rgt'
    cpus 32
    memory '72 GB'
    time '12 h'
=
    input:
    tuple val(meta), path(beds), path(bams)

    output:
    tuple val(meta), path("$differential/*")

    script:
    """
    mkdir -p differential
    rgt-hint differential --organism=hg19 --bc --nc {task.cpus} \\
        --mpbs-files=${beds.join(',')} \\
        --reads-files=.${bams.join(',')} \
        --conditions=${bams.map{ bam -> bam.baseName}.join(',')} \
        --output-location=differential
    """
}
