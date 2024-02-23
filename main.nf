
include { ATAC_LOGFC_BED } from './modules/local/atac_logfc_bed'
include { BEDTOOLS_INTERSECT } from './modules/local/bedtools/intersect'
include { GUNZIP } from './modules/nf-core/gunzip'
include { JOIN_PEAKS_PROMOTERS } from './modules/local/join_peaks_promoters'
include { JOIN_ATAC_RNA } from './modules/local/join_atac_rna'
include { ROWBIND } from './modules/local/rowbind'
include { QUARTO_RENDER } from './modules/local/quarto/render'
include { QUARTONOTEBOOK } from './modules/nf-core/quartonotebook/main'
include { CAT_CAT } from './modules/CCBR/cat/cat/main'

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


    ch_bed_combined = CAT_CAT(ch_atac_bed.map{meta, file -> [ [id: 'combined'], file ]}.groupTuple()).file_out

    BEDTOOLS_INTERSECT(ch_atac_bed.combine(ch_promoters_bed))
    | join(input)
    | JOIN_PEAKS_PROMOTERS
    | JOIN_ATAC_RNA

    ch_atac_rna = JOIN_ATAC_RNA.out.tsv
        .map{ meta, file -> file}
        .collect()
        | ROWBIND
    QUARTONOTEBOOK([[id: 'atac_rna_clusters'], file(params.notebook_atac_rna, checkIfExists: true)],
                   ch_atac_rna.map{file -> [ 'atac_rna_tsv': file.toString() ]},
                   ch_atac_rna,
                   []
    )
}
