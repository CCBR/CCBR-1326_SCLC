process MATRIX_BED {
    """
    Convert a matrix of peak counts to a consensus bed file
    """

    input:
        path(matrix)

    output:
        path("${matrix.baseName}.bed")

    script:
    infile = matrix
    outfile = "${matrix.baseName}.bed"
    template 'matrix_to_bed.R'
}
