
process COUNT_INTERSECT {
    """
    Combine counts from reciprocal bedtools intersect counts
    """
    container 'nciccbr/consensus_peaks:v1.1'

    input:
    tuple val(meta), path(atac_count), path(promoter_count)

    output:
    tuple val(meta), path("counts.${meta.id}.tsv"), emit: tsv

    script:
    outfile = "counts.${meta.id}.tsv"
    template 'count_sets.R'
}
