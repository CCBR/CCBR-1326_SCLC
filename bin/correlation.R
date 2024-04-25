#!/usr/bin/env Rscript
library(dplyr)
library(glue)
library(purrr)
library(readr)
library(stringr)
library(tidyr)

main <- function(peak_gene_infile) {
  dat <- read_tsv(peak_gene_infile)
  metadat <- str_match(peak_gene_infile, "matched_(?<gene>\\w+)_(?<peak>.*)\\.tsv")
  gene_name <- metadat[1, "gene"]
  peak_coord <- metadat[1, "peak"]
  outfile <- str_replace(peak_gene_infile, "\\.tsv$", "_corr.tsv")

  broom::tidy(cor.test(
    dat$rna_count, # escape dollar signs for nextflow template
    dat$atac_count,
    method = "pearson",
    adjust = "fdr"
  )) %>%
    mutate(
      gene_name = gene_name,
      peak_coord = peak_coord
    ) %>%
    write_tsv(outfile)
}

args <- commandArgs(trailingOnly = TRUE)
main(args[1])
