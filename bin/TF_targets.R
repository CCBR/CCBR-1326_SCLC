library(clusterProfiler)
library(dplyr)
library(forcats)
library(ggplot2)
library(glue)
library(here)
library(msigdbr)
library(purrr)
library(readr)
library(rlang)
library(stringr)
library(tibble)
library(tidyr)

eval_parse <- function(x, cd) {
  p <- parse(text = x)
  eval(p, envir = cd)
}

get_gene_list <- function(dat, logfc_col = log2fc_rna, gene_name_col = TSS_gene_name) {
  dat_filt <- dat %>%
    filter(!is.na({{ gene_name_col }}), !is.na({{ logfc_col }})) %>%
    group_by(TSS_gene_name) %>%
    slice_max({{ logfc_col }}, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc({{ logfc_col }}))
  gene_list <- dat_filt %>%
    pull({{ logfc_col }})
  names(gene_list) <- dat_filt %>% pull({{ gene_name_col }})
  return(gene_list)
}

run_gsea <- function(msig_cat_subcat,
                     dat = dat_filt,
                     log2fc_col = log2fc_rna,
                     alpha_level = 0.05, n_top = 15) {
  message(glue("Running GSEA for {msig_cat_subcat} on {as_name(enquo(log2fc_col))}"))
  msig_split <- msig_cat_subcat %>%
    str_split(" ") %>%
    unlist()
  msig_cat <- msig_split[1]
  if (length(msig_split) > 1) {
    msig_sub <- msig_split[2]
    msigdb_set <- msigdbr(species = "Homo sapiens", category = msig_cat, subcategory = msig_sub)
  } else { # no subcategory
    msigdb_set <- msigdbr(species = "Homo sapiens", category = msig_cat)
  }

  term_gene <- msigdb_set %>%
    select(gs_name, gene_symbol)

  # GSEA per cluster
  clusters <- dat %>%
    pull(cluster) %>%
    unique()
  gsea_results <- clusters %>% map(\(cluster_id) {
    gene_list <- dat %>%
      filter(cluster == cluster_id) %>%
      get_gene_list(logfc_col = {{ log2fc_col }})
    return(GSEA(gene_list, TERM2GENE = term_gene, pvalueCutoff = 0.05))
  })
  names(gsea_results) <- clusters

  clusters %>% map(\(cluster_id) {
    try({
      gsea_dotplot <- dotplot(gsea_results[[cluster_id]]) + ggtitle(glue("GSEA {msig_cat_subcat} {cluster_id} {as_name(enquo(log2fc_col))}"))
      ggsave(
        filename = glue("figures/gsea_dotplot_{msig_cat_subcat}_{cluster_id}_{as_name(enquo(log2fc_col))}.png"),
        plot = gsea_dotplot,
        height = 5,
        width = 8
      )
    })
  })
}

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
  select(peak_coord, log2FoldChange, padj, cluster) %>%
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
  select(gene_name, log2FoldChange, padj, cluster) %>%
  rename(
    log2fc_rna = log2FoldChange,
    padj_rna = padj,
    TSS_gene_name = gene_name
  )


dat <- corr_dat %>%
  left_join(rna_dat, relationship = "many-to-many") %>%
  left_join(atac_dat) %>%
  # left_join(motif_dat) %>%
  select(
    TSS_gene_name, peak_coord,
    # TF_list,
    cluster, log2fc_atac, log2fc_rna, padj_atac, padj_rna
  ) %>%
  filter(
    !is.na(log2fc_atac), !is.na(log2fc_rna)
  )

msig_cats <- c("C3 TFT:GTRD", "C3 TFT:TFT_Legacy")
# GSEA on RNA log2fc
msig_cats %>%
  map(\(msig) run_gsea(msig, dat = rna_dat, log2fc_col = log2fc_rna))
# GSEA on ATAC log2fc
msig_cats %>%
  map(\(msig) run_gsea(msig,
    dat = atac_dat %>%
      left_join(motif_dat) %>%
      filter(!is.na(TSS_gene_name), !is.na(log2fc_atac)),
    log2fc_col = log2fc_atac
  ))
