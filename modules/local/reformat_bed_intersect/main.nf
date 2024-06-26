process REFORMAT_BED_INTERSECT_TSS {
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
    tuple val(meta), path("*.reformat_tss.bed")

    script:
    infile = bed
    outfile = "${bed.baseName}.reformat_tss.bed"
    template 'reformat_bed_intersect_TSS.R'

    stub:
    """
    touch ${bed.baseName}.reformat_tss.bed
    """
}
process REFORMAT_BED_INTERSECT_MOTIF {
    """
    Reformat output of bed intersect from peaks bed vs TF motifs bed
    """
    tag { meta.id }
    container "nciccbr/sclc_r-quarto:v0.4.0"
    cpus 1
    memory '16 GB'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.reformat_motif.tsv")

    script:
    infile = bed
    outfile = "${bed.baseName}.reformat_motif.tsv"
    meta_id = meta.id
    template 'reformat_bed_intersect_motif.R'

    stub:
    """
    touch ${bed.baseName}.reformat_motif.tsv
    """
}
