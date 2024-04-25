#!/usr/bin/env Rscript
library(tidyverse)

adjacency <- c(
  "ccbr_tobias/network/c1_vs_c2_c3/c1/adjacency.txt",
  "ccbr_tobias/network/c3_vs_c1_c2/c3/adjacency.txt",
  "ccbr_tobias/network/c2_vs_c1_c3/c2/adjacency.txt"
) %>%
  map(function(f) {
    read_tsv(f) %>%
      mutate(
        cluster_id = str_replace_all(
          f,
          regex("ccbr_tobias/network/.*/(.*)/adjacency.txt"),
          "\\1"
        ),
        outdegree = 1 + str_count(Targets, ",")
      )
  }) %>%
  list_rbind()

write_tsv(adjacency, file = "ccbr_tobias/adjacency_outdegree.tsv")

all_targets <- adjacency %>%
  pull(Targets) %>%
  str_split(",") %>%
  unlist() %>%
  str_remove_all(" ") %>%
  unique() %>%
  sort()

adj_mat <- adjacency %>% select(Source, cluster_id, Targets)


# O_diff

# for tobias we have just one network per group of samples,
# so we don't normalize by number of samples in the group
calc_outdegree_diff <- function(dat, score_col = log2_outdegree_avg) {
  clusters <- dat %>%
    pull(cluster_id) %>%
    unique()
  clusters %>%
    map(function(cluster) {
      z_score_avg_in <- dat %>%
        filter(cluster_id == cluster) %>%
        pull({{ score_col }})
      z_score_avg_out <- dat %>%
        filter(cluster_id != cluster) %>%
        pull({{ score_col }}) %>%
        sum()
      return(tibble(
        cluster_id = cluster,
        score_diff = z_score_avg_in - z_score_avg_out
      ))
    }) %>%
    bind_rows()
}

outdegree_sums <- adjacency %>%
  group_by(cluster_id) %>%
  summarize(sum_outdegree = sum(outdegree))

outdegree_scores <- adjacency %>%
  left_join(cluster_counts) %>%
  left_join(outdegree_sums) %>%
  mutate(log2_outdegree_avg = log2(outdegree / sum_outdegree)) %>%
  group_by(Source) %>%
  nest() %>%
  mutate(O_diff = map(data, calc_outdegree_diff, score_col = log2_outdegree_avg)) %>%
  unnest(O_diff) %>%
  select(-data) %>%
  rename(
    O_diff = score_diff,
    gene_name = Source
  )
