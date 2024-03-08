
include { ATAC_LOGFC_BED           } from './modules/local/atac_logfc_bed'
include { BEDTOOLS_INTERSECT       } from './modules/local/bedtools/intersect'
include { GUNZIP                   } from './modules/nf-core/gunzip'
include { JOIN_PEAKS_PROMOTERS     } from './modules/local/join_peaks_promoters'
include { JOIN_ATAC_RNA            } from './modules/local/join_atac_rna'
include { COUNT_INTERSECT          } from './modules/local/count_intersect'
include { ROWBIND as ROWBIND_COUNT;
          ROWBIND as ROWBIND_LOGFC } from './modules/local/rowbind'
include { CAT_CAT                  } from './modules/CCBR/cat/cat'
include { QUARTONOTEBOOK           } from './modules/nf-core/quartonotebook'
include { MATRIX_BED               } from './modules/local/matrix_bed'
include { CHROMVAR                 } from './modules/local/chromvar'
include { HINT_FOOTPRINTING     } from './modules/local/rgt/hint/footprinting'

workflow {
    ch_consensus_bed = Channel.fromPath(file(params.consensus_peak_matrix, checkIfExists: true)) |
        MATRIX_BED

    ch_cluster_map = Channel.fromPath(file(params.clusters, checkIfExists: true))

    ch_bam = ch_cluster_map
        .splitCsv(header: true, sep: '\t')
        .map{ it ->
            bam = file(it.bam, checkIfExists: true)
            bai = file("${bam}.bai", checkIfExists: true)
            [ [ id: it.sampleName, cluster: it.clusterName ], bam, bai ]
        }

    CHROMVAR(ch_consensus_bed, ch_cluster_map, ch_bam.map{meta, bam, bai -> bam}.collect())

    RGT(ch_consensus_bed, ch_bam)
}

workflow RGT {
    take:
        ch_consensus_bed
        ch_bam

    main:
        ch_rgtdata = Channel.fromPath(file(params.rgtdata)).collect()
        HINT_FOOTPRINTING(ch_bam.combine(ch_consensus_bed), ch_rgtdata)


}

workflow differential {
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

    BEDTOOLS_INTERSECT(ch_atac_bed.combine(ch_promoters_bed)).intersect
    | join(input)
    | JOIN_PEAKS_PROMOTERS
    | JOIN_ATAC_RNA

    COUNT_INTERSECT(BEDTOOLS_INTERSECT.out.counts).tsv
        .map{meta, file -> file}
        .collect()
        .map{ files -> [ [id: 'set_counts'], files ] }
        | ROWBIND_COUNT

    ch_atac_rna = JOIN_ATAC_RNA.out.tsv
        .map{ meta, file -> file}
        .collect()
        .map{ files -> [ [id: 'atac_rna'], files ]}
        | ROWBIND_LOGFC

    qmd_params = ch_atac_rna.combine(ROWBIND_COUNT.out.tsv)
        .map{file1, file2 -> [ 'atac_rna_tsv': file1.toString(), 'set_counts_tsv': file2.toString() ]}
    qmd_inputs = ch_atac_rna.mix(ROWBIND_COUNT.out.tsv).collect()
    QUARTONOTEBOOK([[id: 'atac_rna_clusters'], file(params.notebook_atac_rna, checkIfExists: true)],
                   qmd_params,
                   qmd_inputs,
                   []
    )

}
