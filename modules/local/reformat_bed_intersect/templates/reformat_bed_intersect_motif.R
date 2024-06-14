#!/usr/bin/env Rscript
library(dplyr)
library(readr)
main <- function(infile = "${infile}", outfile = "${outfile}", meta_id = "${meta_id}") {
  read_tsv(infile, col_names = c("peak_chrom", "peak_start", "peak_end", "peak_name", "peak_score", "peak_strand", "TF_chrom", "TF_start", "TF_end", "TF_gene_name", "TF_score", "TF_strand", "overlap")) %>%
    select(peak_chrom, peak_start, peak_end, TF_gene_name, peak_score, peak_strand) %>%
    rename(
      chrom = peak_chrom,
      start = peak_start,
      end = peak_end,
      gene_name = TF_gene_name,
      score = peak_score,
      strand = peak_strand
    ) %>%
    mutate(
      peak_coord = glue("{chr}:{start}-{end}"),
      atac_id_sample = meta_id
    ) %>%
    write_tsv(outfile)
}
main()
