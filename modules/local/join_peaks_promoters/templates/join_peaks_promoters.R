#!/usr/bin/env Rscript
library(dplyr)
library(glue)
library(readr)
library(tidyr)
main <- function(atac_logfc = "${atac_logfc}", atac_promoters = "${atac_promoters}", outfile = "${outfile}") {
  dat_logfc <- read_tsv(atac_logfc) %>% rename(peak_chr_name = `...1`)
  dat_promoters <- read_tsv(atac_promoters,
    col_names = c(
      "peak_chr", "peak_start", "peak_end", "peak_chr_name", "peak_score", "peak_strand",
      "promoter_chr", "promoter_start", "promoter_end", "promoter_chr_name", "promoter_score", "promoter_strand",
      "n_bases_overlap"
    )
  ) %>%
    separate_wider_delim(promoter_chr_name, delim = "|", names = c("promoter_chr_name", "ensemble", "gene_name"))

  dat_promoters %>%
    left_join(dat_logfc, by = "peak_chr_name") %>%
    group_by(peak_chr_name) %>%
    slice_max(n_bases_overlap) %>% # best matching promoter kept for each peak
    select(-c(peak_chr, peak_start, peak_end, peak_score, peak_strand, promoter_chr, promoter_start, promoter_end, promoter_score)) %>%
    write_tsv(outfile)
}
main()
# main(atac_logfc = "data/rank3.pdx.cluster1.alone.logFC_1.pval_0.05.csv",
#      atac_promoters = "output/bedtools_intersect/intersect.rank3.pdx.cluster1.alone.logFC_1.pval_0.05.hg19.promoters.bed",
#      outfile = 'tmp.tsv')
