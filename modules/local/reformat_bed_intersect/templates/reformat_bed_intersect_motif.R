#!/usr/bin/env Rscript
library(dplyr)
library(glue)
library(readr)
main <- function(infile = "${infile}", outfile = "${outfile}", meta_id = "${meta_id}") {
  read_tsv(infile, col_names = c("TF_chrom", "TF_start", "TF_end", "TF_gene_name", "TF_score", "TF_strand", "peak_chrom", "peak_start", "peak_end", "TSS_gene_name", "peak_score", "peak_strand", "overlap")) %>%
    select(peak_chrom, peak_start, peak_end, TF_gene_name, peak_score, peak_strand, TSS_gene_name) %>%
    rename(
      chrom = peak_chrom,
      start = peak_start,
      end = peak_end,
      TF_gene_name = TF_gene_name,
      score = peak_score,
      strand = peak_strand
    ) %>%
    mutate(
      peak_coord = glue("{chrom}:{start}-{end}"),
      atac_sample_id = meta_id
    ) %>%
    distinct() %>%
    write_tsv(outfile)
}
main()
