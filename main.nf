// nf-core
include { GUNZIP                   } from './modules/nf-core/gunzip'
include { QUARTONOTEBOOK           } from './modules/nf-core/quartonotebook'
// CCBR
include { CAT_CAT                  } from './modules/CCBR/cat/cat'
// local
include { ATAC_LOGFC_BED           } from './modules/local/atac_logfc_bed'
include { BEDTOOLS_INTERSECT       } from './modules/local/bedtools/intersect'
include { BEDTOOLS_INTERSECT_FULL as BEDTOOLS_INTERSECT_TSS    } from './modules/local/bedtools/intersect'
include { BEDTOOLS_INTERSECT_FULL as BEDTOOLS_INTERSECT_MOTIF
                                   } from './modules/local/bedtools/intersect'
include { JOIN_PEAKS_PROMOTERS     } from './modules/local/join_peaks_promoters'
include { JOIN_ATAC_RNA            } from './modules/local/join_atac_rna'
include { COUNT_INTERSECT          } from './modules/local/count_intersect'
include { ROWBIND as ROWBIND_COUNT;
          ROWBIND as ROWBIND_LOGFC } from './modules/local/rowbind'
include { MATRIX_BED               } from './modules/local/matrix_bed'
include { CHROMVAR                 } from './modules/local/chromvar/chromvar'
include { CHROMVAR_SUBSET          } from './modules/local/chromvar/chromvar_subset'
include { HINT_FOOTPRINTING        } from './modules/local/rgt/hint/footprinting'
include { MOTIFANALYSIS_MATCHING   } from './modules/local/rgt/motifanalysis'
include { MOTIF2GENE_MAPPING       } from './modules/local/motif2gene_mapping'
include { CLUSTERPROFILER          } from './modules/local/clusterProfiler'
include { EXTRACT_TSS              } from './modules/local/extract_tss'
include { STRIP_TAB                } from './modules/local/strip_tab'
include { FILTER_FOOTPRINTS        } from './modules/local/filter_footprints'
include { REFORMAT_BED_INTERSECT   } from './modules/local/reformat_bed_intersect'
include { NETWORK_OUTDEGREE        } from './modules/local/network_outdegree'
include { CORRELATE_PEAKS_GENES    } from './modules/local/correlate_peaks_genes'

workflow {
    gtf = file(params.gtf, checkIfExists: true)
    pfm = file(params.pfm, checkIfExists: true)
    chrom_sizes = file(params.chrom_sizes, checkIfExists: true)

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

    MOTIF2GENE_MAPPING(gtf, pfm)

    bam_list = ch_bam.map{meta, bam, bai -> bam}.collect()

    // chromvar on all samples
    //CHROMVAR(ch_consensus_bed, ch_cluster_map, bam_list)
    // chromvar on each cluster individually
    //CHROMVAR_SUBSET(ch_consensus_bed.combine(ch_cluster_map).combine(Channel.of('c1', 'c2', 'c3')), bam_list)

    ch_footprints = RGT(ch_consensus_bed, ch_bam).footprints

    ch_tss_bed = EXTRACT_TSS(gtf, chrom_sizes).bed

    ch_peaks_tss = BEDTOOLS_INTERSECT_TSS(ch_consensus_bed.map{
                bed -> [[id: 'consensus'], bed]
            }.combine(ch_tss_bed)
        ).intersect
        | REFORMAT_BED_INTERSECT
        | map{ meta, bed -> bed }

    Channel.fromPath([file(params.metadata, checkIfExists: true),
                      file(params.rna_counts_norm, checkIfExists: true),
                      file(params.atac_counts_norm, checkIfExists: true)
                     ])
        .collect()
        .combine(ch_peaks_tss)
        | CORRELATE_PEAKS_GENES

    ch_tf_genes = BEDTOOLS_INTERSECT_MOTIF(ch_footprints.combine(ch_peaks_tss)).intersect
        | FILTER_FOOTPRINTS
    ch_tf_genes | NETWORK_OUTDEGREE
    ß
}

workflow RGT {
    take:
        ch_consensus_bed
        ch_bam

    main:

        ch_rgtdata = Channel.fromPath(file(params.rgtdata)).collect()
        /*
        // HINT_FOOTPRINTING keeps rerunning even though prior runs were successful
        */
        // HINT_FOOTPRINTING(ch_bam.combine(ch_consensus_bed), ch_rgtdata)
        ch_hint_beds = Channel.fromPath("output/hint_footprinting/*_footprints/*.bed") |
            map{ bed ->
                bed_id = bed.baseName
                [ [id: bed_id], bed ]
            }
        /*
        // keeps rerunning even though prior runs were successful...

        MOTIFANALYSIS_MATCHING(
            ch_hint_beds, //HINT_FOOTPRINTING.out.bed,
            ch_rgtdata
        )
        */
        ch_footprints = Channel.fromPath("output/motifanalysis_matching/*_motifs/*.bed") |
            map{ bed ->
                bed_id = bed.baseName.replace("_mpbs", "")
                [ [id: bed_id], bed ]
            } |
            STRIP_TAB

    emit:
        footprints = STRIP_TAB.out

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

    ch_atac_rna | CLUSTERPROFILER

    qmd_params = ch_atac_rna.combine(ROWBIND_COUNT.out.tsv)
        .map{file1, file2 -> [ 'atac_rna_tsv': file1.toString(), 'set_counts_tsv': file2.toString() ]}
    qmd_inputs = ch_atac_rna.mix(ROWBIND_COUNT.out.tsv).collect()
    QUARTONOTEBOOK([[id: 'atac_rna_clusters'], file(params.notebook_atac_rna, checkIfExists: true)],
                   qmd_params,
                   qmd_inputs,
                   []
    )

}
