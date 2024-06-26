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
# parallel(
#   gene_name = 'MIEF1',
#   input_files = 'matched_MIEF1_chr22:39893805-39895071.tsv,matched_MIEF1_chr22:39900282-39900340.tsv,matched_MIEF1_chr22:39897088-39898779.tsv,matched_MIEF1_chr22:39892405-39892628.tsv' %>%
#     str_split(',') %>% unlist() %>% paste0('output_2/split_pairs/matches/', .) %>% str_c(collapse = ','),
#   num_cores = 8
#   )
