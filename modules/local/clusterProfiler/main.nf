process CLUSTERPROFILER {
    cpus 16
    memory '16 GB'
    time '24 h'

    container "nciccbr/sclc_r-quarto:v0.2.1"

    input:
        path(tsv)

    output:
        path('clusterProfiler.RData'),  emit: rda
        path('enrich_ORA_results.png'), emit: png

    script:
    output_rda = 'clusterProfiler.RData'
    output_png = 'enrich_ORA_results.png'
    template 'run_clusterProfiler.R'
}
