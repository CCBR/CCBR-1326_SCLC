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
    tuple val(meta), path("*.filt.tsv")

    script:
    outfile = "${infile.baseName}.filt.tsv"
    template 'filter_footprints.R'

    stub:
    """
    touch ${infile.baseName}.filt.bed
    """
}
