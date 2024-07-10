process CHROMVAR {
    cpus 16
    memory '16 GB'
    time '24 h'

    container "nciccbr/sclc_r-quarto:v0.4.1"

    input:
        path(consensus_bed)
        path(cluster_map)
        path(bams)

    output:
        path('chromVAR.results*')

    script:
    bam_filenames = bams.join(',')
    output_tsv = 'chromVAR.results.tsv'
    output_rda = 'chromVAR.results.RData'
    output_png = 'chromVAR.results.variability.png'
    template 'chromvar.R'
}
