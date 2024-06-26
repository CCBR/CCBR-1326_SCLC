#!/usr/bin/env Rscript
options(error = rlang::entrace)
library(dplyr)
library(furrr)
library(future)
library(glue)
library(purrr)
library(readr)
library(stringr)
library(tidyr)

main <- function(infile = "output_2/peaks_genes_motifs/peaks_genes_motifs.tsv",
                 num_cores = 8) {
  message(glue("using {num_cores} cores for parallel processing"))
  plan(multicore, workers = num_cores)
  dat <- read_tsv(infile)

  dir.create("matches", showWarnings = FALSE)
  nest_dat <- dat %>%
    group_by(peak_coord, TSS_gene_name) %>%
    nest() %>%
    head(n = 20) # TODO delete later -- this is for debugging with a smaller dataset
  nest_dat %>%
    future_pmap(\(TSS_gene_name, peak_coord, data) {
      filename <- glue("matches/matched_{TSS_gene_name}_{peak_coord}.tsv")
      data %>%
        mutate(TSS_gene_name = TSS_gene_name, peak_coord = peak_coord) %>%
        write_tsv(filename)
    })
}

args <- commandArgs(trailingOnly = TRUE)
main(infile = args[1], num_cores = as.integer(args[2]))
