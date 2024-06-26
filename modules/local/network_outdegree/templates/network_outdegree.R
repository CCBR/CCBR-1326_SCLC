#!/usr/bin/env Rscript
library(dplyr)
library(readr)

main <- function(samples_infile = "${matched_samples}",
                 corr_infile = "${correlations}",
                 outfile = "${outfile}") {
  read_tsv(samples_infile) %>%
    right_join(read_tsv(corr_infile)) %>%
    separate_longer_delim(TF_list, ",") %>%
    rename(TF_name = TF_list) %>%
    select(TF_name, TSS_gene_name, atac_sample_id) %>%
    distinct() %>%
    group_by(TF_name) %>%
    summarize(outdegree = n()) %>%
    write_tsv(outfile)
}
main()
