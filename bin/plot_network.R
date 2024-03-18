library(dplyr)
library(readr)

library(ggraph)
library(ggiraph)
library(tidygraph)

df <- read_tsv("ccbr_tobias/network/c2_vs_c1_c3/c2/edges.txt") %>%
  rename(source = Sites_3, target = Origin_0)

nodes <- tibble(node = c(df %>% pull(source), df %>% pull(target)) %>%
  unique() %>%
  sort()) %>%
  mutate(rownum = row_number())

edges <- df %>%
  select(source, target) %>%
  distinct() %>%
  left_join(nodes, by = c("source" = "node")) %>%
  rename(from = rownum) %>%
  left_join(nodes, by = c("target" = "node")) %>%
  rename(to = rownum) %>%
  select(from, to)

graph <- as_tbl_graph(edges, nodes = nodes)

p <- ggraph(graph) +
  geom_edge_link() +
  geom_node_point()
girafe(ggobj = p)
