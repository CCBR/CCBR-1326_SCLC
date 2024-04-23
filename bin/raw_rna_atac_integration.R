library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(tidyverse)

# Retrieve gene exon lengths
# https://support.bioconductor.org/p/p132346/#p132372
ex <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene") %>% IRanges::reduce()
exlen <- relist(width(unlist(ex)), ex)
exlens <- sapply(exlen, sum)

# Map entrez gene IDs to HGNC gene symbols
# https://www.biostars.org/p/69647/#69648
gene_symbols <- annotate::getSYMBOL(names(exlens), data = "org.Hs.eg")
genes_df <- tibble(hgnc_symbol = gene_symbols, entrez_id = names(gene_symbols)) %>%
  right_join(tibble(entrez_id = names(exlens), gene_length = exlens),
    by = "entrez_id"
  ) %>%
  mutate(gene_len_kb = gene_length / 1000)

# Calculate TMM-normalized RPKM with edgeR
raw_rna_counts <- read_tsv("data/RNA.rawcount.subsetted.protein.coding_sample_name_modified.tsv")
gene_lengths_df <- genes_df %>%
  right_join(raw_rna_counts, by = c("hgnc_symbol" = "Gene_Id")) %>%
  dplyr::select(hgnc_symbol, gene_length) %>%
  dplyr::rename(Length = gene_length)
rna_edger <- edgeR::DGEList(
  counts = raw_rna_counts %>% dplyr::select(-Gene_Id),
  genes = gene_lengths_df
) %>%
  edgeR::calcNormFactors(method = "TMM")

rna_counts_norm <- bind_cols(gene_lengths_df, edgeR::rpkm(rna_edger))

write_csv(rna_counts_norm, "data/rna_counts_normalized.csv")
