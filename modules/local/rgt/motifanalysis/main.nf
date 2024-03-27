process MOTIFANALYSIS_MATCHING {
    //container "nciccbr/sclc_rgt:v0.1.0"
    module 'rgt'
    cpus 32
    memory '2 GB'
    time '12 h'

    input:
    tuple val(meta), path(bed)
    path(rgtdata)

    output:
    tuple val(meta), path("${meta.id}_motifs/*_mpbs.bed"), emit: bed
    tuple val(meta), path("${meta.id}_motifs/*"),                   emit: motif_dir

    script:
    outfile_rgt = "${meta.id}_motifs/${meta.id}_mpbs.bed"
    """
    mkdir ${meta.id}_motifs/
    rgt-motifanalysis matching \\
        --organism=hg19 \\
        --output-location=${meta.id}_motifs/ \\
        --input-files ${bed}
    """
}
