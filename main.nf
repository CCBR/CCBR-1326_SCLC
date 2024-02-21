
include { ATAC_LOGFC_BED } from './modules/local/atac_logfc_bed'

workflow {
    input = Channel.fromPath(file(params.datasheet, checkIfExists: true))
        .splitCsv(header: true)
        .map {
            [ [id: "c${it.cluster}"], file(it.atac, checkIfExists: true), file(it.rna, checkIfExists: true) ]
        }
    input
        .map{ meta, atac, rna -> [ meta, atac ] }
        | ATAC_LOGFC_BED
        | view
}
