process MOTIF2GENE_MAPPING {
    """
    Map TFs to the genes that encode them
    """
    cpus 1

    container 'nciccbr/ccbr_ubuntu_base_20.04:v6'

    input:
    path(gtf)
    path(pfm)

    output:
    path('*.motif2gene_mapping.txt')

    script:
    output_txt = "${pfm.baseName}_${gtf.baseName}.motif2gene_mapping.txt"
    template 'motif2gene_mapping.py'
}
