process CHROMVAR_SUBSET {
    cpus 16
    memory '16 GB'
    time '24 h'

    container "nciccbr/sclc_r-quarto:v0.2.1"

    input:
        path(consensus_bed)
        path(cluster_map)
        path(bams)
        val(cluster_id)

    output:
        path("chromVAR.${cluster_id}*")

    script:
    bam_filenames = bams.join(',')
    output_tsv = "chromVAR.${cluster_id}.tsv"
    output_rda = "chromVAR.${cluster_id}.RData"
    output_png = "chromVAR.${cluster_id}.variability.png"
    template 'chromvar.R'
}
