#!/usr/bin/env Rscript
library(dplyr)
library(readr)
library(tidyr)

main <- function(input_tsv = "output_2/correlations/gene_peak_corr.tsv",
                 output_bed = "tmp.bed",
                 pvalue_thresh = 0.01) {
  read_tsv(input_tsv) %>%
    filter(estimate > 0, p.value < pvalue_thresh) %>%
    separate_wider_delim(peak_coord, delim = ":", names = c("chr", "pos")) %>%
    separate_wider_delim(pos, delim = "-", names = c("start", "end")) %>%
    select(chr, start, end, TSS_gene_name) %>%
    filter(!is.na(chr), !is.na(start), !is.na(end), !is.na(TSS_gene_name))
  mutate(score = ".", strand = ".") %>%
    arrange(chr, start, end) %>%
    write_tsv(output_bed, col_names = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
main(args[1], args[2])
