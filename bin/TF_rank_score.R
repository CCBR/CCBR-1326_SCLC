#!/usr/bin/env Rscript
# code based on supplemental methods of 10.1126/science.abe1505
library(tidyverse)

cluster_membership <- read_tsv("assets/cluster_membership.tsv") %>%
  rename(
    cluster_id = clusterName,
    sample_id = sampleName
  )

cluster_counts <- cluster_membership %>%
  dplyr::count(cluster_id)

cluster_membership %<>%
  left_join(cluster_counts, by = "cluster_id") %>%
  select(-bam)


# A_diff (Accessibility difference) was calculated as the average z-score difference between one group relative to others
chromvar_dat <- read_tsv("output/chromvar/chromVAR.results.tsv")

calc_zscore_diff <- function(dat, score_col = z_score) {
  clusters <- dat %>%
    pull(cluster_id) %>%
    unique()
  dat_sum <- dat %>%
    group_by(cluster_id) %>%
    summarise(z_score_sum = sum({{ score_col }}))
  clusters %>%
    map(function(cluster) {
      z_score_avg_in <- dat_sum %>%
        filter(cluster_id == cluster) %>%
        pull(z_score_sum) / cluster_counts %>%
          filter(cluster_id == cluster) %>%
          pull(n)
      z_score_avg_out <- dat_sum %>%
        filter(cluster_id != cluster) %>%
        pull(z_score_sum) %>%
        sum() / cluster_counts %>%
          filter(cluster_id == cluster) %>%
          pull(n) %>%
          sum()
      return(tibble(cluster_id = cluster, score_diff = z_score_avg_in - z_score_avg_out))
    }) %>%
    bind_rows()
}

atac_scores <- chromvar_dat %>%
  pivot_longer(ends_with(".bam"),
    names_to = "sample_id", values_to = "z_score"
  ) %>%
  mutate(sample_id = str_replace_all(sample_id, ".sorted.dedup.bam", "")) %>%
  left_join(cluster_membership) %>%
  select(name, sample_id, z_score, cluster_id, n) %>%
  group_by(name) %>%
  nest() %>%
  mutate(A_diff = map(data, calc_zscore_diff, score_col = z_score)) %>%
  unnest(A_diff) %>%
  select(-data) %>%
  rename(A_diff = score_diff)


# E_diff (Expression difference) as the third criterion,
# based on the assumption that TFs with higher relative expression are more important in that subtype of samples
diff_dat <- read_tsv("output/rowbind_logfc/concat.atac_rna.tsv")

rna_scores <- diff_dat %>%
  filter(pvalue_rna < 0.05) %>%
  select(gene_name, log2FoldChange_rna, pvalue_rna, cluster_id) %>%
  # rename_with( ~str_remove(.x, '_rna'), .cols = ends_with('_rna')) %>%
  mutate(
    minus_log10_pvalue = -log10(pvalue_rna),
    E_diff = minus_log10_pvalue * log2FoldChange_rna / abs(log2FoldChange_rna)
  ) %>%
  select(gene_name, cluster_id, E_diff)

# O_diff
adjacency <- read_tsv("ccbr_tobias/adjacency_outdegree.tsv")

adjacency
