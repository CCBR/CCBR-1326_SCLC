#!/usr/bin/env Rscript
library(tidyverse)
main <- function(peak_gene_infile = "${tsv}") {
    dat <- read_tsv(peak_gene_infile)
    metadat <- str_match(peak_gene_infile, "(?<gene>\\\\w+)_(?<peak>.*)\\\\.tsv")
    gene_name <- metadat[1,'gene']
    peak_coord <- metadat[1,'peak']
    outfile <- str_replace(peak_gene_infile, '\\\\.tsv\$', '_corr.tsv')

    broom::tidy(cor.test(
        .x\$rna_count, # escape dollar signs for nextflow template
        .x\$atac_count,
        method = "spearman",
        adjust = "fdr"
    )) %>%
        mutate(gene_name = gene_name,
               peak_coord = peak_coord) %>%
        write_tsv(outfile)
}
main()
