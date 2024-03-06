process MATRIX_BED {
    """
    Convert a matrix of peak counts to a consensus bed file
    """
    container "nciccbr/sclc_r-quarto:v0.2.1"

    input:
        path(matrix)

    output:
        path("${matrix.baseName}.bed")

    script:
    infile = matrix
    outfile = "${matrix.baseName}.bed"
    template 'matrix_to_bed.R'
}
