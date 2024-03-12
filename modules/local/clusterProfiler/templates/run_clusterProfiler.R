library(clusterProfiler)
library(dplyr)
library(msigdbr)
library(readr)

tsv_filename <- "output/rowbind_logfc/concat.atac_rna.tsv" # "${tsv}"

# C2 CP database from msigdbr
# TODO try C5 C6 C7
term_gene <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP") %>%
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

results <- dat %>%
  pivot_longer(c(log2FoldChange_atac, log2FoldChange_rna),
    names_to = "log2fc_source", values_to = "log2FoldChange"
  ) %>%
  mutate(log2fc_source = str_remove(log2fc_source, "log2FoldChange_")) %>%
  group_by(cluster_id, log2fc_source) %>%
  nest() %>%
  mutate(
    gene_list = map(data, get_gene_list),
    enrich = map(gene_list, function(gl) {
      enricher(names(unlist(gl)), TERM2GENE = term_gene)
    }),
    gsea = map(gene_list, function(gl) {
      GSEA(unlist(gl), TERM2GENE = term_gene)
    })
  )

saveRDS(results, "clusterProfiler_results.Rds")

genes_c2_rna <- dat %>%
  filter(cluster_id == "c2") %>%
  get_gene_list(logfc_col = log2FoldChange_rna)
enrich_c2cp <- enricher(names(genes_c2_rna), TERM2GENE = term_gene)
gsea_c2cp <- GSEA(genes_c2_rna, TERM2GENE = term_gene)
