list.files("data", pattern = "*.edit.tsv", full.names = TRUE) %>%
  map(\(x) read_tsv(x) %>% mutate(filename = x)) %>%
  bind_rows() %>%
  mutate(cluster_id = str_extract(filename, "cluster\\d", ) %>%
    str_remove("luster")) %>%
  select(gene_name, cluster_id, log2FoldChange, padj) %>%
  pivot_wider(
    names_from = cluster_id,
    values_from = c(log2FoldChange, padj)
  ) %>%
  write_tsv("data/rna_log2fc_for_IPA.tsv")
