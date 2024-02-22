
include { ATAC_LOGFC_BED } from './modules/local/atac_logfc_bed'
include { BEDTOOLS_INTERSECT } from './modules/local/bedtools/intersect'
include { GUNZIP } from './modules/nf-core/gunzip'

workflow {
    input = Channel.fromPath(file(params.datasheet, checkIfExists: true))
        .splitCsv(header: true)
        .map {
            [ [id: "c${it.cluster}"], file(it.atac, checkIfExists: true), file(it.rna, checkIfExists: true) ]
        }
    ch_atac_bed = input
        .map{ meta, atac, rna -> [ meta, atac ] }
        | ATAC_LOGFC_BED

    
    Channel.fromPath(file(params.promoters, checkIfExists: true))
        .map{ bed -> 
            [ [id: 'promoters'], bed ]
        }
        | GUNZIP
    ch_promoters_bed = GUNZIP.out.gunzip.map{ meta, bed -> bed }
    
    BEDTOOLS_INTERSECT(ch_atac_bed.combine(ch_promoters_bed))
}
   
