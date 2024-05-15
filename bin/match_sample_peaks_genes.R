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
           atac_counts_infile = "data/raw_tmm_fpkm_batch_corrected_PDX_only.csv",
           peak_gene_infile = "output/reformat_bed_intersect/intersect.raw_tmm_fpkm_batch_corrected_PDX_only.gencode.v19.annotation.TSS_padded.reformat.bed",
           matched_outfile = "matched_genes_peaks.tsv") {
    # integrate RNA-seq with ATAC-seq
    pdx_metadat <- read_tsv(pdx_meta_infile) %>%
      rename(sample_id_verbose = sample_id) %>%
      mutate(
        sample_id = str_extract(sample_id_verbose,
          pattern = "Sample_\\d{1,2}_([a-zA-Z0-9]{6})_.*",
          group = 1
        )
      ) %>%
      select(sample_id, sample_id_verbose, ATAC_Id, rank3) %>%
      rename_with(~ tolower(glue("pdx_{.x}")))
    metadat <- readxl::read_excel(metadata_infile) %>%
      mutate(sample_id = as.character(Sample_ID)) %>%
      full_join(pdx_metadat, by = c("sample_id" = "pdx_sample_id"))

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
        names_to = "atac_id_sample",
        values_to = "atac_count"
      ) %>%
      mutate(atac_sample_id = str_extract(atac_id_sample,
        pattern = "([a-zA-Z0-9]+)_S.*",
        group = 1
      )) %>%
      right_join(peak_gene_tss, relationship = "many-to-many")

    atac_ids <- atac_counts_norm_long %>%
      select(atac_id_sample, atac_sample_id) %>%
      distinct()

    rna_counts_norm_long <- rna_counts_norm %>%
      rename(gene_name = hgnc_symbol) %>%
      right_join(peak_gene_tss) %>%
      pivot_longer(-c(peak_coord, gene_name),
        names_to = "rna_id_sample",
        values_to = "rna_count"
      ) %>%
      mutate(rna_sample_id = str_extract(rna_id_sample,
        pattern = "Sample_\\d{1,2}_([a-zA-Z0-9]{6})_.*",
        group = 1
      ))

    rna_ids <- rna_counts_norm_long %>%
      select(rna_id_sample, rna_sample_id) %>%
      distinct()

    metadat_join <- metadat %>%
      rename(
        rna_id_meta = `Matching Bulk RNA ID`,
        atac_id_meta = ATAC_ID
      ) %>%
      mutate(
        rna_id_meta = case_when(
          rna_id_meta == "N" &
            !is.na(`Parth Found_RNA_ID`) ~ `Parth Found_RNA_ID`,
          rna_id_meta == "N" ~ NA_character_,
          TRUE ~ rna_id_meta
        )
      ) %>%
      # filter(!is.na(rna_id_meta)) %>%
      select(pdx_sample_id_verbose, sample_id, atac_id_meta, rna_id_meta, pdx_rank3) %>%
      separate_longer_delim(rna_id_meta, "/") %>%
      mutate(atac_id = str_extract(atac_id_meta,
        pattern = "([a-zA-Z0-9]+)",
        group = 1
      )) %>%
      full_join(atac_ids,
        by = c("atac_id" = "atac_sample_id")
      ) %>%
      full_join(rna_ids,
        by = c("pdx_sample_id_verbose" = "rna_id_sample")
      ) %>%
      select(-rna_sample_id) %>%
      rename(
        rna_sample_id = pdx_sample_id_verbose,
        rna_id = sample_id,
        atac_sample_id = atac_id_sample
      ) %>%
      filter(atac_id != "2705020")
    samples_mapped <- metadat_join %>%
      filter(!is.na(atac_id), !is.na(rna_id))

    # free up memory before `counts_join`
    rm(
      list = c(
        "atac_counts_norm",
        "metadat",
        "peak_gene_tss",
        "rna_counts_norm"
      )
    )
    counts_join <- full_join(
      rna_counts_norm_long %>%
        rename(rna_id = rna_sample_id) %>%
        right_join(samples_mapped, relationship = "many-to-many"),
      atac_counts_norm_long %>%
        rename(atac_id = atac_sample_id) %>%
        right_join(samples_mapped, relationship = "many-to-many"),
      relationship = "many-to-many"
    )
    # free up memory after `counts_join`
    rm(
      list = c(
        "atac_counts_norm_long",
        "rna_counts_norm_long",
        "samples_mapped"
      )
    )
    counts_sum <- counts_join %>%
      filter(!is.na(atac_count), !is.na(rna_count)) %>%
      select(-rna_sample_id, ends_with("_orig")) %>%
      group_by(gene_name, peak_coord) %>%
      summarize(n = n()) %>%
      filter(n > 10)

    counts_join %>%
      filter(!is.na(atac_count), !is.na(rna_count)) %>%
      right_join(counts_sum %>% select(-n), # only keep peak-gene pairs that are in at least 10 samples
        by = c("gene_name", "peak_coord")
      ) %>%
      write_tsv(matched_outfile)

    # counts_grp <- counts_join %>%
    #   filter(!is.na(atac_count), !is.na(rna_count)) %>%
    #   right_join(counts_sum %>% select(-n), # only keep peak-gene pairs that are in at least 10 samples
    #     by = c("gene_name", "peak_coord")
    #   ) %>%
    #   select(
    #     -rna_sample_id,
    #     rna_id_sample, rna_id, rna_count, atac_id, atac_count
    #   ) %>%
    #   group_by(gene_name, peak_coord) %>%
    #   nest()

    # head(counts_grp)

    ## error: creates too many output files for nextflow to copy
    ### /usr/bin/ls: Argument list too long
    # l <- counts_grp %>%
    #   pmap(\(gene_name, peak_coord, data) {
    #     filename <- glue("matched_{gene_name}_{peak_coord}.tsv")
    #     write_tsv(data, filename)
    #   })
  }


args <- commandArgs(trailingOnly = TRUE)
main(
  metadata_infile = args[1],
  pdx_meta_infile = args[2],
  rna_counts_infile = args[3],
  atac_counts_infile = args[4],
  peak_gene_infile = args[5]
)
