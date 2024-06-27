#!/usr/bin/env Rscript
# code based on supplemental methods of 10.1126/science.abe1505
library(tidyverse)
library(rlang)
library(patchwork)
theme_sovacool <- ggplot2::theme_bw() +
  ggplot2::theme(
    legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(1, "lines"),
    plot.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
    strip.background = element_blank(),
    strip.placement = "outside"
  )
theme_set(theme_sovacool)

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
      z_score_avg_in <- (dat_sum %>%
        filter(cluster_id == cluster) %>%
        pull(z_score_sum)) / (cluster_counts %>%
        filter(cluster_id == cluster) %>%
        pull(n))
      z_score_avg_out <- (dat_sum %>%
        filter(cluster_id != cluster) %>%
        pull(z_score_sum) %>%
        sum()) / (cluster_counts %>%
        filter(cluster_id != cluster) %>%
        pull(n) %>%
        sum())
      return(tibble(
        cluster_id = cluster,
        score_diff = z_score_avg_in - z_score_avg_out
      ))
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
  rename(
    A_diff = score_diff,
    gene_name = name
  )


# E_diff (Expression difference)
# based on the assumption that TFs with higher relative expression are more important in that subtype of samples
diff_dat <- read_tsv("output/rowbind_logfc/concat.atac_rna.tsv")

expr_scores <- diff_dat %>%
  select(gene_name, log2FoldChange_rna, pvalue_rna, cluster_id) %>%
  mutate(
    minus_log10_pvalue = -log10(pvalue_rna),
    E_diff = minus_log10_pvalue * log2FoldChange_rna / abs(log2FoldChange_rna)
  ) %>%
  select(gene_name, cluster_id, E_diff)


# O_diff (Out-degree difference)
# We assume if a TF is important in one subtype of samples,
# it would regulate more genes in that subtype relative to other subtypes,
# as demonstrated with a higher O_diff

outdegree_dat <- read_tsv("output_2/network_outdegree/TF_outdegree.tsv") %>%
  left_join(cluster_membership)

outdegree_sums <- outdegree_dat %>%
  group_by(sample_id) %>%
  nest() %>%
  mutate(sum_outdegree = map(data, ~ sum(.x$outdegree))) %>%
  select(-data) %>%
  unnest(sum_outdegree) %>%
  ungroup()
outdegree_dat_norm <- outdegree_dat %>%
  left_join(cluster_counts) %>%
  left_join(outdegree_sums) %>%
  mutate(log2_outdegree_norm = log2(outdegree / sum_outdegree))
outdegree_scores <- outdegree_dat_norm %>%
  group_by(TF_name) %>%
  nest() %>%
  mutate(O_diff = map(data, calc_zscore_diff, score_col = log2_outdegree_norm)) %>%
  unnest(O_diff) %>%
  select(-data) %>%
  rename(
    O_diff = score_diff,
    gene_name = TF_name
  )

# overall rank
rank_col <- function(dat, value_col = A_diff, group_col = cluster_id) {
  dat %>%
    group_by({{ group_col }}) %>%
    arrange(desc({{ value_col }})) %>%
    mutate("rank_{{value_col}}" := row_number())
}

dat_joined <- atac_scores %>%
  full_join(expr_scores, by = c("gene_name", "cluster_id")) %>%
  full_join(outdegree_scores, by = c("gene_name", "cluster_id"))

tf_rank_dat <- dat_joined %>%
  filter(if_all(ends_with("_diff"), ~ !is.na(.))) %>%
  rank_col(value_col = A_diff) %>%
  rank_col(value_col = E_diff) %>%
  rank_col(value_col = O_diff) %>%
  mutate(rank_sum = rank_A_diff + rank_E_diff + rank_O_diff) %>%
  group_by(cluster_id) %>%
  arrange(rank_sum) %>%
  mutate(TF_rank = row_number()) %>%
  arrange(TF_rank)


plot_heatmap_row <- function(dat, value_column = A_diff,
                             scale_fill = scale_fill_viridis_c) {
  dat %>%
    ggplot(aes(
      x = gene_name, y = 1, # y = diff,
      fill = {{ value_column }}
    )) +
    geom_tile() +
    facet_wrap(~cluster_id, nrow = 1, scales = "free") +
    labs(
      y = quo_name(enquo(value_column)),
      x = ""
    ) +
    theme(
      axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
    )
}
top_tfs <- tf_rank_dat %>%
  filter(TF_rank <= 25) %>%
  mutate(
    outdegree = O_diff,
    expression = E_diff,
    accessibility = A_diff
  )

write_csv(top_tfs, "output_2/top_TFs_HINT.csv")

min_max <- list(
  min = ~ min(.x, na.rm = TRUE),
  max = ~ max(.x, na.rm = TRUE)
)
tf_rank_dat %>%
  ungroup() %>%
  summarize(across(ends_with("_diff"), min_max))
# plot TFs patchwork
(plot_heatmap_row(top_tfs, outdegree) +
  colorspace::scale_fill_continuous_diverging(
    palette = "Blue-Red",
  ) +
  theme(axis.text.x = element_blank())
) /
  (plot_heatmap_row(top_tfs, expression) +
    colorspace::scale_fill_continuous_diverging(
      palette = "Green-Orange",
    ) +
    theme(
      axis.text.x = element_blank(),
      strip.text = element_blank()
    )
  ) /
  (plot_heatmap_row(top_tfs, accessibility) +
    colorspace::scale_fill_continuous_diverging(
      palette = "Purple-Green",
    ) +
    theme(strip.text = element_blank()))

# debugging distribution of outdegrees
outdegree_dat %>%
  ggplot(aes(outdegree, fill = cluster_id)) +
  geom_histogram(alpha = 0.5, position = position_identity(), bins = 20)

outdegree_dat_norm %>%
  ggplot(aes(log2_outdegree_norm, fill = cluster_id)) +
  geom_histogram(alpha = 0.5, position = position_identity(), bins = 20)

outdegree_dat_norm %>%
  select(sample_id, sum_outdegree, cluster_id) %>%
  distinct() %>%
  ggplot(aes(sum_outdegree, fill = cluster_id)) +
  geom_histogram(alpha = 0.5, position = position_identity())

outdegree_scores %>%
  ggplot(aes(O_diff, fill = cluster_id)) +
  geom_histogram(alpha = 0.5, position = position_identity(), bins = 20)

dat <- outdegree_dat %>%
  left_join(cluster_counts) %>%
  left_join(outdegree_sums) %>%
  mutate(log2_outdegree_norm = log2(outdegree / sum_outdegree)) %>%
  filter(TF_name == "ARNT2")
dat_sum <- dat %>%
  group_by(cluster_id) %>%
  summarise(z_score_sum = sum(log2_outdegree_norm))
cluster <- "c3"
z_score_avg_in <- (dat_sum %>%
  filter(cluster_id == cluster) %>%
  pull(z_score_sum)) / (cluster_counts %>%
  filter(cluster_id == cluster) %>%
  pull(n))
z_score_avg_out <- (dat_sum %>%
  filter(cluster_id != cluster) %>%
  pull(z_score_sum) %>%
  sum()) / (cluster_counts %>%
  filter(cluster_id != cluster) %>%
  pull(n) %>%
  sum())
score_diff <- z_score_avg_in - z_score_avg_out

tf_rank_dat %>%
  pivot_longer(c(A_diff, E_diff, O_diff)) %>%
  ggplot(aes(value, fill = cluster_id)) +
  geom_histogram(alpha = 0.5, position = position_identity(), bins = 20) +
  facet_wrap(~name, scales = "free")

write_tsv(tf_rank_dat, file = "output_2/tf_rank_HINT.tsv")
