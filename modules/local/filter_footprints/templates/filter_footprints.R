#!/usr/bin/env Rscript
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

main <- function(infile = "${infile}", outfile = "${outfile}") {
  read_tsv(infile, col_names = c("TF_chrom", "TF_start", "TF_end", "motif_name", "TF_score", "TF_strand", "peak_chrom", "peak_start", "peak_end", "peak_gene_name", "peak_score", "peak_strand", "n_bases_overlap")) %>%
    mutate(peak_width = TF_end - TF_start) %>%
    filter(TF_score > 10, peak_width >= 10, peak_width <= 20) %>%
    select(-peak_width) %>%
    mutate(TF_name = str_remove(motif_name, "MA\\\\d+.\\\\d.") %>% str_remove("\\\\(.*\\\\)")) %>%
    separate_longer_delim(TF_name, "::") %>%
    write_tsv(outfile, col_names = TRUE)
}
main()
