process REFORMAT_BED_INTERSECT {
    """
    Reformat output of bed intersect from consensus bed vs TSS bed
    to keep the gene name and strand info from the TSS bed file
    """
    tag { meta.id }
    container "nciccbr/sclc_r-quarto:v0.3.3"
    cpus 1
    memory '16 GB'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.reformat.bed")

    script:
    infile = bed
    outfile = "${bed.baseName}.reformat.bed"
    template 'reformat_bed_intersect.R'

    stub:
    """
    touch ${bed.baseName}.reformat.bed
    """
}
