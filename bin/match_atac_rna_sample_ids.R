#!/usr/bin/env Rscript
options(error = rlang::entrace)
library(dplyr)
library(glue)
library(purrr)
library(readr)
library(stringr)
library(tidyr)

main <-
  function(metadata_infile = "assets/metadata_for_kelley.xlsx",
           pdx_meta_infile = "assets/pdx_rank3_metadata.tsv",
           rna_counts_infile = "data/rna_counts_normalized.csv",
           atac_counts_infile = "data/raw_tmm_fpkm_batch_corrected_PDX_only.csv") {
    # integrate RNA-seq with ATAC-seq
    pdx_metadat <- read_tsv(pdx_meta_infile) %>%
      rename(sample_id_verbose = sample_id) %>%
      mutate(
        sample_id = str_extract(sample_id_verbose, pattern = "Sample_\\d{1,2}_([a-zA-Z0-9]{6})_.*", group = 1)
      ) %>%
      select(sample_id, sample_id_verbose, ATAC_Id, rank3) %>%
      rename_with(~ tolower(glue("pdx_{.x}")))
    metadat <- readxl::read_excel(metadata_infile) %>%
      mutate(sample_id = as.character(Sample_ID)) %>%
      full_join(pdx_metadat, by = c("sample_id" = "pdx_sample_id"))

    rna_counts_norm <- read_csv(rna_counts_infile)
    atac_counts_norm <-
      read_csv(atac_counts_infile)

    atac_counts_norm_long <- atac_counts_norm %>%
      rename(peak_coord = Coordinate) %>%
      pivot_longer(-peak_coord, names_to = "atac_id_sample", values_to = "atac_count") %>%
      mutate(atac_sample_id = str_extract(atac_id_sample, pattern = "([a-zA-Z0-9]+)_S.*", group = 1))

    atac_ids <- atac_counts_norm_long %>%
      select(atac_id_sample, atac_sample_id) %>%
      distinct()

    rna_counts_norm_long <- rna_counts_norm %>%
      rename(gene_name = hgnc_symbol) %>%
      pivot_longer(-gene_name, names_to = "rna_id_sample", values_to = "rna_count") %>%
      mutate(
        rna_sample_id = str_extract(rna_id_sample, pattern = "Sample_\\d{1,2}_([a-zA-Z0-9]{6})_.*", group = 1)
      )

    rna_ids <- rna_counts_norm_long %>%
      select(rna_id_sample, rna_sample_id) %>%
      distinct()

    metadat_join <- metadat %>%
      rename(rna_id_meta = `Matching Bulk RNA ID`, atac_id_meta = ATAC_ID) %>%
      mutate(
        rna_id_meta = case_when(
          rna_id_meta == "N" &
            !is.na(`Parth Found_RNA_ID`) ~ `Parth Found_RNA_ID`,
          rna_id_meta == "N" ~ NA_character_,
          TRUE ~ rna_id_meta
        )
      ) %>%
      # filter(!is.na(rna_id_meta)) %>%
      select(
        pdx_sample_id_verbose,
        sample_id,
        atac_id_meta,
        rna_id_meta,
        pdx_rank3
      ) %>%
      separate_longer_delim(rna_id_meta, "/") %>%
      mutate(atac_id = str_extract(atac_id_meta, pattern = "([a-zA-Z0-9]+)", group = 1)) %>%
      full_join(atac_ids, by = c("atac_id" = "atac_sample_id")) %>%
      full_join(rna_ids, by = c("pdx_sample_id_verbose" = "rna_id_sample")) %>%
      select(-rna_sample_id) %>%
      rename(
        rna_sample_id = pdx_sample_id_verbose,
        rna_id = sample_id,
        atac_sample_id = atac_id_sample
      ) %>%
      filter(atac_id != "2705020")

    samples_mapped <- metadat_join %>%
      filter(!is.na(atac_id), !is.na(rna_id))
    samples_mapped %>%
      write_tsv("matched_sample_ids.tsv")
    rna_counts_norm_long %>%
      rename(rna_id = rna_sample_id) %>%
      write_tsv("rna_counts_norm_long.tsv")
    atac_counts_norm_long %>%
      rename(atac_id = atac_sample_id) %>%
      write_tsv("atac_counts_norm_long.tsv")
  }


args <- commandArgs(trailingOnly = TRUE)
main(
  metadata_infile = args[1],
  pdx_meta_infile = args[2],
  rna_counts_infile = args[3],
  atac_counts_infile = args[4]
)
