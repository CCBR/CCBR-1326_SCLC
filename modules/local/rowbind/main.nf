
process ROWBIND {
    """
    Bind rows of multiple tsv files
    """
    container 'nciccbr/consensus_peaks:v1.1'

    input:
    path(input_list)

    output:
    path("concat.atac_rna.tsv"), emit: tsv

    script:
    infiles_list = input_list.join(',')
    outfile = "concat.atac_rna.tsv"
    template 'rowbind.R'
}
