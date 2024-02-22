#!/usr/bin/env Rscript
library(dplyr)
library(glue)
library(readr)
library(stringr)
library(tidyr)
main <- function(promoters_logfc = "${promoters_logfc}", rna_logfc = "${rna_logfc}", outfile = "${outfile}") {
  dat_rna <- read_csv(rna_logfc) %>% rename(gene_name = `...1`)
  dat_promoters <- read_tsv(promoters_logfc)
  genes_atac_rna <- dat_promoters %>%
    inner_join(dat_rna, by = "gene_name") %>%
    rename_with(~ str_replace(., ".x", "_atac"), ends_with(".x")) %>%
    rename_with(~ str_replace(., ".y", "_rna"), ends_with(".y"))
  message(glue("RNAseq:\t{nrow(dat_rna)}",
    "ATACseq:\t{nrow(dat_promoters)}",
    "joined:\t{nrow(genes_atac_rna)}",
    .sep = "\n"
  ))
  genes_atac_rna %>% write_tsv(outfile)
}
main()
# main(promoters_logfc = "output/join_peaks_promoters/intersect.rank3.pdx.cluster1.alone.logFC_1.pval_0.05.hg19.promoters.tsv",
#      rna_logfc = "data/cluster1_vs_other_DGE.csv",
#     outfile = 'tmp.tsv')
