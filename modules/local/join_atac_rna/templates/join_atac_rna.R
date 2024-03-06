#!/usr/bin/env Rscript
library(dplyr)
library(glue)
library(readr)
library(stringr)
library(tidyr)
main <- function(promoters_logfc = "${promoters_logfc}",
                 rna_logfc = "${rna_logfc}",
                 outfile = "${outfile}",
                 logfile = "${logfile}",
                 cluster_id = "${meta.id}") {
  dat_rna <- read_tsv(rna_logfc) %>% rename(gene_name = `...1`)
  dat_promoters <- read_tsv(promoters_logfc)
  genes_atac_rna <- dat_promoters %>%
    inner_join(dat_rna, by = "gene_name") %>%
    rename_with(~ str_replace(., ".x", "_atac"), ends_with(".x")) %>%
    rename_with(~ str_replace(., ".y", "_rna"), ends_with(".y")) %>%
    group_by(gene_name) %>%
    slice_max(n_bases_overlap) %>% # pick one peak per gene based on number of bases that overlap
    ungroup()
  genes_atac_rna %>%
    mutate(cluster_id = cluster_id) %>%
    write_tsv(outfile)

  log_msg <- glue("RNAseq:\t{nrow(dat_rna)}",
    "ATACseq:\t{nrow(dat_promoters)}",
    "joined:\t{nrow(genes_atac_rna)}",
    .sep = "\n"
  )
  write_lines(log_msg, logfile)
}
main()
# main(promoters_logfc = "output/join_peaks_promoters/intersect.rank3.pdx.cluster1.DEseq2.hg19.promoters.atac-promoters.tsv",
#      rna_logfc = "data/rank3.pdx.cluster1_RNA.DEseq2.tsv",
#     outfile = 'tmp.tsv', logfile = 'tmp.log')
