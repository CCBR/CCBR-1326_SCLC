#!/usr/bin/env Rscript
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

main <- function(infile = "${infile}", outfile = "${outfile}") {
  read_tsv(infile, col_names = c("TF_chrom", "TF_start", "TF_end", "motif_name", "TF_score", "TF_strand")) %>%
    mutate(peak_width = TF_end - TF_start) %>%
    filter(TF_score > 10, peak_width >= 10, peak_width <= 20) %>%
    select(-peak_width) %>%
    mutate(TF_name = str_remove(motif_name, "MA\\\\d+.\\\\d.") %>% str_remove("\\\\(.*\\\\)")) %>%
    separate_longer_delim(TF_name, "::") %>%
    select(TF_chrom, TF_start, TF_end, TF_name, TF_score, TF_strand) %>%
    write_tsv(outfile, col_names = FALSE)
}
main()
