#!/usr/bin/env Rscript

library(dplyr)
library(glue)
library(readr)
library(tidyr)

main <- function(cluster_id = "${meta.id}", atac_count = "${atac_count}", promoter_count = "${promoter_count}", outfile = "${outfile}") {
  read_tsv(atac_count, col_names = c("chr", "start", "end", "chr_start_end", "score", "strand", "count")) %>%
    select(chr_start_end, count) %>%
    mutate(count_type = "peak_per_gene") %>%
    bind_rows(read_tsv(promoter_count, col_names = c("chr", "start", "end", "chr_start_end", "score", "strand", "count")) %>%
      select(chr_start_end, count) %>%
      mutate(count_type = "gene_per_peak")) %>%
    mutate(cluster_id = cluster_id) %>%
    write_tsv(outfile)
}

main()
# main(
#   cluster_id = "c1",
#   atac_count = "output/bedtools_intersect/count.atac-promoter.rank3.pdx.cluster1.alone.logFC_1.pval_0.05.hg19.promoters.bed",
#   promoter_count = "output/bedtools_intersect/count.promoter-atac.hg19.promoters.rank3.pdx.cluster1.alone.logFC_1.pval_0.05.bed",
#   outfile = "temp.tsv"
# )
