#!/usr/bin/env Rscript
library(dplyr)
library(readr)

main <- function(infile = "${infile}", outfile = "${outfile}") {
  tf_gene_links <- read_tsv(infile)
  tf_gene_links %>%
    select(motif_name, TF_name, peak_gene_name) %>%
    distinct() %>%
    group_by(TF_name) %>%
    summarize(outdegree = n()) %>%
    mutate(sample_id = "${meta.id}") %>%
    write_tsv(outfile)
}
main()
