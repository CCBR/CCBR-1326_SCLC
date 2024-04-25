#!/usr/bin/env Rscript
library(dplyr)
library(glue)
library(purrr)
library(readr)
library(stringr)
library(tidyr)

main <-
  function(metadata_infile = "assets/matched_atac_RNAdata.csv",
           rna_counts_infile = "data/rna_counts_normalized.csv",
           atac_counts_infile = "data/raw_tmm_fpkm_batch_corrected_PDX_only.csv",
           peak_gene_infile = "output/reformat_bed_intersect/intersect.raw_tmm_fpkm_batch_corrected_PDX_only.gencode.v19.annotation.TSS_padded.reformat.bed") {
    # integrate RNA-seq with ATAC-seq
    metadat <- read_csv(metadata_infile)
    rna_counts_norm <- read_csv(rna_counts_infile)
    atac_counts_norm <-
      read_csv(atac_counts_infile)
    peak_gene_tss <-
      read_tsv(
        peak_gene_infile,
        col_names = c("chr", "start", "end", "gene_name", "score", "strand")
      ) %>%
      mutate(peak_coord = glue("{chr}:{start}-{end}")) %>%
      select(gene_name, peak_coord)

    atac_counts_norm_long <- atac_counts_norm %>%
      rename(peak_coord = Coordinate) %>%
      pivot_longer(-peak_coord,
        names_to = "atac_sample_id_orig",
        values_to = "atac_count"
      ) %>%
      mutate(
        atac_sample_id = str_extract(atac_sample_id_orig,
          pattern = "([a-zA-Z0-9]+)_S.*",
          group = 1
        )
      ) %>%
      right_join(peak_gene_tss, relationship = "many-to-many")

    atac_ids <- atac_counts_norm_long %>%
      select(atac_sample_id_orig, atac_sample_id) %>%
      distinct()

    rna_counts_norm_long <- rna_counts_norm %>%
      rename(gene_name = hgnc_symbol) %>%
      right_join(peak_gene_tss) %>%
      pivot_longer(-c(peak_coord, gene_name),
        names_to = "rna_sample_id",
        values_to = "rna_count"
      )

    rna_ids <- rna_counts_norm_long %>%
      select(rna_sample_id) %>%
      distinct()

    metadat_join <- metadat %>%
      filter(`Matching Bulk RNA (Y/N)` == "Y") %>%
      rename(
        rna_id = `Matching Bulk RNA ID`,
        atac_id_orig = ATAC_ID
      ) %>%
      select(atac_id_orig, rna_id) %>%
      separate_longer_delim(rna_id, "/") %>%
      mutate(rna_id = str_remove(rna_id, "Sample_")) %>%
      mutate(atac_id = str_extract(atac_id_orig,
        pattern = "([a-zA-Z0-9]+)",
        group = 1
      )) %>%
      full_join(
        atac_ids,
        by = c("atac_id" = "atac_sample_id")
      ) %>%
      full_join(
        rna_ids,
        by = c("rna_id" = "rna_sample_id")
      )
    samples_mapped <- metadat_join %>%
      filter(!is.na(atac_id), !is.na(rna_id)) %>%
      rename(
        rna_sample_id = rna_id,
        atac_sample_id = atac_id
      )

    counts_join <- full_join(
      rna_counts_norm_long %>%
        right_join(samples_mapped, relationship = "many-to-many"),
      atac_counts_norm_long %>%
        right_join(samples_mapped, relationship = "many-to-many"),
      relationship = "many-to-many"
    )
    # free up memory
    rm(
      list = c(
        "atac_counts_norm",
        "atac_counts_norm_long",
        "metadat",
        "metadat_join",
        "peak_gene_tss",
        "rna_counts_norm",
        "rna_counts_norm_long",
        "samples_mapped"
      )
    )
    counts_sum <- counts_join %>%
      filter(!is.na(atac_count), !is.na(rna_count)) %>%
      select(-rna_sample_id) %>%
      group_by(gene_name, peak_coord) %>%
      summarize(n = n()) %>%
      filter(n > 10)

    counts_join %>%
      filter(!is.na(atac_count) | !is.na(rna_count)) %>%
      right_join(counts_sum %>% select(-n), # only keep peak-gene pairs that are in at least 10 samples
        by = c("gene_name", "peak_coord")
      ) %>%
      select(atac_sample_id, rna_sample_id) %>%
      distinct()


    counts_grp <- counts_join %>%
      filter(!is.na(atac_count), !is.na(rna_count)) %>%
      right_join(counts_sum %>% select(-n), # only keep peak-gene pairs that are in at least 10 samples
        by = c("gene_name", "peak_coord")
      ) %>%
      select(-rna_sample_id) %>%
      group_by(gene_name, peak_coord) %>%
      nest()

    counts_grp %>%
      pmap(\(gene_name, peak_coord, data) {
        filename <- glue("matched_{gene_name}_{peak_coord}.tsv")
        write_tsv(data, filename)
      })
  }


args <- commandArgs(trailingOnly = TRUE)
main(
  metadata_infile = args[1],
  rna_counts_infile = args[2],
  atac_counts_infile = args[3],
  peak_gene_infile = args[4]
)
