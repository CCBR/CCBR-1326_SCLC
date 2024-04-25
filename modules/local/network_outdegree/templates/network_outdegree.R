#!/usr/bin/env Rscript
library(dplyr)
library(readr)

main <- function(infile = "${infile}", outfile = "${outfile}") {
  tf_gene_links <- read_tsv(infile, col_names = c(
    "TF_chrom", "TF_start", "TF_end", "TF_name", "TF_score", "TF_strand",
    "peak_chrom", "peak_start", "peak_end", "peak_gene_name", "peak_score", "peak_strand",
    "n_bases_overlap"
  ))
  tf_gene_links %>%
    select(TF_name, peak_gene_name) %>%
    distinct() %>%
    group_by(TF_name) %>%
    summarize(outdegree = n()) %>%
    mutate(sample_id = "${meta.id}") %>%
    write_tsv(outfile)
}
main()
