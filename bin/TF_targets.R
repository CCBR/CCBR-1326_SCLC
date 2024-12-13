library(tidyverse)
library(here)
alpha_level <- 0.05
log2fc_thresh <- 1

corr_dat <- read_tsv(here("output_2", "filter_corr_tsv", "gene_peak_corr_filt.tsv"))
motif_dat <- read_tsv(here("output_2", "peaks_genes_motifs", "peaks_genes_motifs.tsv"))

atac_dat <- "data/rank3.pdx.cluster1.DEseq2.tsv data/rank3.pdx.cluster2.DEseq2.tsv data/rank3.pdx.cluster3.DEseq2.tsv" %>%
  str_split(" ") %>%
  unlist() %>%
  map(\(filename) {
    read_tsv(filename) %>%
      mutate(
        atac_filename = filename,
        cluster = str_extract(filename, "cluster[1-3]")
      ) %>%
      rename(peak_coord = `...1`)
  }) %>%
  bind_rows() %>%
  select(peak_coord, log2FoldChange, pvalue, padj, cluster) %>%
  rename(
    log2fc_atac = log2FoldChange,
    padj_atac = padj
  )

rna_dat <- "data/rank3.pdx.cluster1_RNA.DEseq2.edit.tsv data/rank3.pdx.cluster3_RNA.DEseq2.edit.tsv data/rank3.pdx.cluster2_RNA.DEseq2.edit.tsv" %>%
  str_split(" ") %>%
  unlist() %>%
  map(\(filename) {
    read_tsv(filename) %>%
      mutate(
        rna_filename = filename,
        cluster = str_extract(filename, "cluster[1-3]")
      )
  }) %>%
  bind_rows() %>%
  select(gene_name, log2FoldChange, pvalue, padj, cluster) %>%
  rename(
    log2fc_rna = log2FoldChange,
    padj_rna = padj,
    TSS_gene_name = gene_name
  )


dat <- corr_dat %>%
  left_join(motif_dat) %>%
  left_join(atac_dat) %>%
  left_join(rna_dat) %>%
  select(TSS_gene_name, peak_coord, TF_list, cluster, atac_sample_id, rna_sample_id, atac_count, rna_count, log2fc_atac, log2fc_rna, padj_atac, padj_rna)
dat_filt <- dat %>%
  filter(
    !is.na(log2fc_atac), !is.na(log2fc_rna),
    padj_atac < alpha_level, abs(log2fc_atac) >= log2fc_thresh,
    padj_rna < alpha_level, abs(log2fc_rna) >= log2fc_thresh
  )

dat_filt %>%
  select(-atac_sample_id, rna_sample_id, atac_count, rna_count) %>%
  separate_longer_delim(cols = TF_list, delim = ",") %>%
  rename(TF_name = TF_list)


dat_filt %>%
  select(TSS_gene_name, peak_coord, cluster) %>%
  unique()
dat %>%
  pull(TSS_gene_name) %>%
  unique() %>%
  length()
rna_dat %>%
  pull(TSS_gene_name) %>%
  unique() %>%
  length()
corr_dat %>%
  left_join(motif_dat) %>%
  left_join(atac_dat) %>%
  pull(TSS_gene_name) %>%
  unique() %>%
  length()
