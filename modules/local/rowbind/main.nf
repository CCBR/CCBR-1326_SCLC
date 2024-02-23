
process ROWBIND {
    """
    Bind rows of multiple tsv files
    """
    container 'nciccbr/consensus_peaks:v1.1'

    input:
    tuple val(meta), path(input_list)

    output:
    path("concat.*.tsv"), emit: tsv

    script:
    infiles_list = input_list.join(',')
    outfile = "concat.${meta.id}.tsv"
    template 'rowbind.R'
}
