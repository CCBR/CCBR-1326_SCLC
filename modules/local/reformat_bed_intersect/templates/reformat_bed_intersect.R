#!/usr/bin/env Rscript
library(dplyr)
library(readr)
main <- function(infile = "${infile}", outfile = "${outfile}") {
  read_tsv(infile, col_names = c("peak_chrom", "peak_start", "peak_end", "peak_name", "peak_score", "peak_strand", "TSS_chrom", "TSS_start", "TSS_end", "TSS_gene_name", "TSS_score", "TSS_strand", "overlap")) %>%
    select(peak_chrom, peak_start, peak_end, TSS_gene_name, peak_score, TSS_strand) %>%
    write_tsv(outfile, col_names = FALSE)
}
main()
