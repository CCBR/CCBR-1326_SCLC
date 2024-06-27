process NETWORK_OUTDEGREE {
    """
    Calculate network outdegree with TF motifs from HINT-ATAC linked with TSS regions
    """
    container "nciccbr/sclc_r-quarto:v0.3.3"
    cpus 1
    memory '32 GB'

    input:
    tuple path(matched_samples), path(correlations)

    output:
    path("TF_outdegree.tsv")

    script:
    outfile = "TF_outdegree.tsv"
    template 'network_outdegree.R'

    stub:
    """
    touch ${infile.baseName}.outdegree.tsv
    """
}
