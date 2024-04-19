library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(tidyverse)

# Retrieve gene exon lengths
# https://support.bioconductor.org/p/p132346/#p132372
ex <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene") %>% reduce()
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

# Calculate FPKM
# https://www.biostars.org/p/312185/#312188
raw_rna_counts <- read_tsv("data/RNA.rawcount.subsetted.protein.coding_sample_name_modified.tsv") %>%
  pivot_longer(!matches("Gene_Id"),
    names_to = "sample_id", values_to = "count"
  )
sample_cpm <- raw_rna_counts %>%
  group_by(sample_id) %>%
  summarize(cpm = sum(count) / 10^6)
raw_rna_counts %>%
  full_join(sample_cpm, by = "sample_id") %>%
  left_join(genes_df, by = c("Gene_Id" = "hgnc_symbol")) %>%
  mutate(fpkm = count / cpm / gene_len_kb)
