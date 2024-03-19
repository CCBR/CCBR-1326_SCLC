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
          regex("ccbr_tobias/network/c1_vs_c2_c3/(.*)/adjacency.txt"),
          "\\1"
        ),
        outdegree = 1 + str_count(Targets, ",")
      )
  }) %>%
  list_rbind()
