process FILTER_FOOTPRINTS {
    """
    Filter TF footprints from RGT
    """
    tag { meta.id }
    container "nciccbr/sclc_r-quarto:v0.3.3"
    cpus 1
    memory '32 GB'

    input:
    tuple val(meta), path(infile)

    output:
    tuple val(meta), path("*.filt.bed")

    script:
    outfile = "${infile.baseName}.filt.bed"
    template 'filter_footprints.R'

    stub:
    """
    touch ${infile.baseName}.filt.bed
    """
}
