library(clusterProfiler)
library(dplyr)
library(msigdbr)
library(readr)
library(rlang)
library(tibble)
library(ggplot2)

# nextflow variables
tsv_filename <- "${tsv}" # "output/rowbind_logfc/concat.atac_rna.tsv"
output_rda <- "${output_rda}" # "clusterProfiler.RData"
output_png <- "${output_png}" # 'enrich_ORA_results.png'

# C5 = Gene Ontology
msigdb_set <- msigdbr(species = "Homo sapiens", category = "C5")
term_gene <- msigdb_set %>%
  select(gs_name, gene_symbol)

alpha_level <- 0.05
log2fc_thresh <- 0 # setting to zero removes filtering by log2FC

dat <- read_tsv(tsv_filename) %>%
  filter(
    pvalue_rna < alpha_level,
    pvalue_atac < alpha_level,
    abs(log2FoldChange_atac) >= log2fc_thresh,
    abs(log2FoldChange_rna) >= log2fc_thresh,
  )

get_gene_list <- function(dat, logfc_col = log2FoldChange) {
  dat_filt <- dat %>%
    filter(!is.na(gene_name)) %>%
    arrange(desc({{ logfc_col }}))
  gene_list <- dat_filt %>% pull({{ logfc_col }})
  names(gene_list) <- dat_filt %>% pull(gene_name)
  return(gene_list)
}

# ORA & GSEA per cluster

dat_long <- dat %>%
  pivot_longer(c(log2FoldChange_atac, log2FoldChange_rna),
    names_to = "log2fc_source", values_to = "log2FoldChange"
  ) %>%
  mutate(log2fc_source = str_remove(log2fc_source, "log2FoldChange_"))

gsea_results <- dat_long %>%
  group_by(cluster_id, log2fc_source) %>%
  nest() %>%
  mutate(
    gene_list = map(data, get_gene_list),
    gsea = map(gene_list, function(gl) {
      GSEA(unlist(gl), TERM2GENE = term_gene)
    })
  )


# evaluate GeneRatio fractions which are character vectors
# https://community.rstudio.com/t/evaluale-strings-as-expression-in-mutate/164655/2
eval_parse <- function(x, cd) {
  p <- parse(text = x)
  eval(p, envir = cd)
}

enrich_results <- dat_long %>%
  filter(log2fc_source == "atac") %>% # enricher just uses list of genes, not log2 values
  group_by(cluster_id) %>%
  nest() %>%
  mutate(
    gene_list = map(data, get_gene_list),
    enrich = map(gene_list, function(gl) {
      enricher(names(unlist(gl)), TERM2GENE = term_gene) %>%
        as_tibble()
    })
  ) %>%
  select(cluster_id, enrich) %>%
  unnest(cols = enrich) %>%
  mutate(category = str_split_i(ID, "_", 1)) %>%
  rowwise() %>%
  mutate(gene_ratio = eval_parse(GeneRatio, pick(everything()))) %>%
  ungroup()

n_top <- 20
top_enrich <- enrich_results %>%
  group_by(cluster_id, ID) %>%
  summarize(max_gr = max(gene_ratio)) %>%
  slice_max(order_by = max_gr, n = n_top)
enrich_plot <- enrich_results %>%
  filter(ID %in% (top_enrich %>% pull(ID))) %>%
  mutate(name = fct_reorder(ID, .x = gene_ratio, .fun = max)) %>%
  ggplot(aes(gene_ratio, name, size = Count, color = cluster_id)) +
  geom_point(alpha = 0.9) + # , position = position_dodge(width = 0.3)) +
  geom_hline(
    yintercept = seq(1.5, length(enrich_results %>% pull(ID) %>% unique()), 1),
    lwd = 0.5, colour = "grey92"
  ) +
  scale_color_viridis_d() +
  # facet_wrap(~ cluster_id, nrow = 3, scales = 'free') +
  theme_bw() +
  theme(panel.grid.major.y = element_blank())
ggsave(filename = output_png, plot = enrich_plot)


save.image(output_rda)
