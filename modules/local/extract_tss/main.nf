process EXTRACT_TSS {
    """
    Extract gene TSSs and write to BED file with ± 0.5Mb of padding
    """
    cpus 1

    container 'nciccbr/ccbr_ubuntu_base_20.04:v6'

    input:
    path(gtf)
    path(chrom_sizes)

    output:
    path('*TSS_padded.bed'), emit: bed

    script:
    bed = "${gtf.baseName}.TSS_padded.bed"
    """
    extract_gene_TSS_bed.py ${gtf} ${chrom_sizes} ${bed}
    """
}
