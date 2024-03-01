#!/usr/bin/env Rscript
library(dplyr)
library(readr)
library(tidyr)
main <- function(infile = "${infile}", outfile = "${outfile}") {
  read_csv(infile) %>%
    rename(chrom_start_end = Coordinate) %>%
    separate_wider_delim(chrom_start_end, delim = ":", names = c("chrom", "start_end"), cols_remove = FALSE) %>%
    separate_wider_delim(start_end, delim = "-", names = c("start", "end")) %>%
    select(chrom, start, end, chrom_start_end) %>%
    mutate(
      score = ".",
      strand = "."
    ) %>%
    write_tsv(outfile, col_names = FALSE)
}
main()
# main(infile = "data/raw_tmm_fpkm_batch_corrected_PDX_only.csv",
#      outfile = "tmp.bed")
