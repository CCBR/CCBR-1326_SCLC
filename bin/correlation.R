#!/usr/bin/env Rscript
library(dplyr)
library(furrr)
library(future)
library(glue)
library(purrr)
library(readr)
library(stringr)
library(tidyr)

correlate <- function(peak_gene_infile) {
  dat <- read_tsv(peak_gene_infile)
  metadat <- str_match(peak_gene_infile, "matched_(?<gene>\\w+)_(?<peak>.*)\\.tsv")
  TSS_gene_name <- metadat[1, "gene"]
  peak_coord <- metadat[1, "peak"]

  corr_result <- broom::tidy(cor.test(
    dat$rna_count,
    dat$atac_count,
    method = "pearson",
    adjust = "fdr"
  )) %>%
    mutate(
      TSS_gene_name = TSS_gene_name,
      peak_coord = peak_coord
    )
  return(corr_result)
}

parallel <- function(gene_name, input_files, num_cores = 12) {
  message(glue("using {num_cores} cores for parallel processing"))
  plan(multicore, workers = num_cores)
  tsv_files <- input_files %>%
    str_split(",") %>%
    unlist()
  tsv_files %>%
    furrr::future_map(correlate) %>%
    bind_rows() %>%
    write_tsv(glue("{gene_name}_corr.tsv"))
}

args <- commandArgs(trailingOnly = TRUE)
parallel(
  gene_name = args[1],
  input_files = args[2],
  num_cores = as.integer(args[3])
)
