process NETWORK_OUTDEGREE {
    """
    Calculate network outdegree with TF motifs from HINT-ATAC linked with TSS regions
    """
    tag { meta.id }
    container "nciccbr/sclc_r-quarto:v0.3.3"
    cpus 1
    memory '32 GB'

    input:
    tuple val(meta), path(infile)

    output:
    tuple val(meta), path("${infile.baseName}.outdegree.tsv")

    script:
    outfile = "${infile.baseName}.outdegree.tsv"
    template 'network_outdegree.R'

    stub:
    """
    touch ${infile.baseName}.outdegree.tsv
    """
}
