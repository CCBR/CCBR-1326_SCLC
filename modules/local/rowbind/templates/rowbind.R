#!/usr/bin/env Rscript
library(dplyr)
library(stringr)
library(purrr)
library(readr)

main <- function(infiles_list = "${infiles_list}", outfile = "${outfile}") {
  str_split(infiles_list, ",") %>%
    purrr::map(read_tsv) %>%
    purrr::list_rbind() %>%
    write_tsv(outfile)
}
main()
