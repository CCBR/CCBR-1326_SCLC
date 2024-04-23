#!/usr/bin/env Rscript
library(tidyverse)
library(glue)

# integrate RNA-seq with ATAC-seq
metadat <- read_csv("assets/matched_atac_RNA_metadata.csv")
rna_counts_norm <- read_csv("data/rna_counts_normalized.csv")
atac_counts_norm <- read_csv("data/raw_tmm_fpkm_batch_corrected_PDX_only.csv")
peak_gene_tss <- read_tsv("output/reformat_bed_intersect/intersect.raw_tmm_fpkm_batch_corrected_PDX_only.gencode.v19.annotation.TSS_padded.reformat.bed",
  col_names = c("chr", "start", "end", "gene_name", "score", "strand")
) %>%
  mutate(peak_coord = glue("{chr}:{start}-{end}")) %>%
  select(gene_name, peak_coord)

atac_counts_norm_long <- atac_counts_norm %>%
  rename(peak_coord = Coordinate) %>%
  pivot_longer(-peak_coord,
    names_to = "atac_sample_id",
    values_to = "atac_count"
  ) %>%
  mutate(atac_sample_id = str_extract(atac_sample_id,
    pattern = "([a-zA-Z0-9]+)[a-zA-Z0-9_]+",
    group = 1
  )) %>%
  right_join(peak_gene_tss, relationship = "many-to-many")


rna_counts_norm_long <- rna_counts_norm %>%
  rename(gene_name = hgnc_symbol) %>%
  right_join(peak_gene_tss) %>%
  pivot_longer(-c(peak_coord, gene_name),
    names_to = "rna_sample_id",
    values_to = "rna_count"
  )

metadat_join <- metadat %>%
  filter(`Matching Bulk RNA (Y/N)` == "Y") %>%
  rename(RNA_ID = `Matching Bulk RNA ID`) %>%
  select(ATAC_ID, RNA_ID) %>%
  separate_longer_delim(RNA_ID, "/") %>%
  mutate(RNA_ID = str_remove(RNA_ID, "Sample_")) %>%
  mutate(ATAC_ID = str_extract(ATAC_ID,
    pattern = "([a-zA-Z0-9]+)[a-zA-Z0-9_]+",
    group = 1
  )) %>%
  full_join(atac_counts_norm_long %>% select(atac_sample_id) %>% distinct(),
    by = c("ATAC_ID" = "atac_sample_id")
  ) %>%
  full_join(rna_counts_norm_long %>% select(rna_sample_id) %>% distinct(),
    by = c("RNA_ID" = "rna_sample_id")
  )
samples_mapped <- metadat_join %>%
  filter(!is.na(ATAC_ID), !is.na(RNA_ID)) %>%
  rename(
    rna_sample_id = RNA_ID,
    atac_sample_id = ATAC_ID
  )

counts_join <- full_join(
  rna_counts_norm_long %>%
    right_join(samples_mapped, relationship = "many-to-many"),
  atac_counts_norm_long %>%
    right_join(samples_mapped, relationship = "many-to-many"),
  relationship = "many-to-many"
)
counts_corr <- counts_join %>%
  filter(!is.na(atac_count), !is.na(rna_count)) %>%
  select(-rna_sample_id) %>%
  group_by(gene_name, peak_coord) %>%
  nest() %>%
  mutate(
    test = map(data, ~ cor.test(.x$rna_count, .x$atac_count,
      method = "spearman", adjust = "fdr"
    )),
    tidied = map(test, tidy)
  ) %>%
  unnest(cols = tidied) %>%
  select(-data, -test)

write_csv(counts_corr, "output/peaks_genes_corr.csv")
