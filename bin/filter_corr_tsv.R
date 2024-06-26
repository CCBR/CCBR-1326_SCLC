#!/usr/bin/env Rscript
library(dplyr)
library(readr)
library(tidyr)

main <- function(input_tsv = "output_2/correlations/gene_peak_corr.tsv",
                 output_tsv = "tmp.tsv",
                 pvalue_thresh = 0.01) {
  read_tsv(input_tsv) %>%
    filter(estimate > 0, p.value < pvalue_thresh) %>%
    arrange(desc(estimate)) %>%
    write_tsv(output_tsv)
}

args <- commandArgs(trailingOnly = TRUE)
main(args[1], args[2])
