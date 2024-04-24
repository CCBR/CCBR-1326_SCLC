#!/usr/bin/env Rscript
library(dplyr)
library(glue)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(future)
library(furrr)

main <-
  function(metadata_infile = "assets/matched_atac_RNA_metadata.csv",
           rna_counts_infile = "data/rna_counts_normalized.csv",
           atac_counts_infile = "data/raw_tmm_fpkm_batch_corrected_PDX_only.csv",
           peak_gene_infile = "output/reformat_bed_intersect/intersect.raw_tmm_fpkm_batch_corrected_PDX_only.gencode.v19.annotation.TSS_padded.reformat.bed",
           peak_gene_outfile = "output/peaks_genes_corr.csv",
           ncores = 16,
           pvalue_thresh = 0.01) {
    n_workers <- min(ncores, availableCores())
    print(glue("Using {n_workers} workers"))
    plan(multisession, workers = 8)

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
        names_to = "atac_sample_id",
        values_to = "atac_count"
      ) %>%
      mutate(
        atac_sample_id = str_extract(atac_sample_id,
          pattern = "([a-zA-Z0-9]+)[a-zA-Z0-9_]+",
          group = 1
        )
      ) %>%
      right_join(peak_gene_tss, relationship = "many-to-many")


    rna_counts_norm_long <- rna_counts_norm %>%
      rename(gene_name = hgnc_symbol) %>%
      right_join(peak_gene_tss) %>%
      pivot_longer(-c(peak_coord, gene_name),
        names_to = "rna_sample_id",
        values_to = "rna_count"
      )

    metadat_join <- metadat %>%
      filter(`Matching Bulk RNA (Y/N)` == "Y") %>%
      rename(RNA_ID = `Matching Bulk RNA ID`) %>%
      select(ATAC_ID, RNA_ID) %>%
      separate_longer_delim(RNA_ID, "/") %>%
      mutate(RNA_ID = str_remove(RNA_ID, "Sample_")) %>%
      mutate(ATAC_ID = str_extract(ATAC_ID,
        pattern = "([a-zA-Z0-9]+)[a-zA-Z0-9_]+",
        group = 1
      )) %>%
      full_join(
        atac_counts_norm_long %>% select(atac_sample_id) %>% distinct(),
        by = c("ATAC_ID" = "atac_sample_id")
      ) %>%
      full_join(
        rna_counts_norm_long %>% select(rna_sample_id) %>% distinct(),
        by = c("RNA_ID" = "rna_sample_id")
      )
    samples_mapped <- metadat_join %>%
      filter(!is.na(ATAC_ID), !is.na(RNA_ID)) %>%
      rename(
        rna_sample_id = RNA_ID,
        atac_sample_id = ATAC_ID
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

    counts_corr <- counts_join %>%
      filter(!is.na(atac_count), !is.na(rna_count)) %>%
      right_join(counts_sum %>% select(-n), # only keep peak-gene pairs that are in at least 10 samples
        by = c("gene_name", "peak_coord")
      ) %>%
      select(-rna_sample_id) %>%
      group_by(gene_name, peak_coord) %>%
      nest() %>%
      mutate(test = future_map(data, ~ broom::tidy(
        cor.test(
          .x\$rna_count, # escape dollar signs for nextflow template
          .x\$atac_count,
          method = "spearman",
          adjust = "fdr"
        )
      ))) %>%
      unnest(cols = test) %>%
      select(-data) %>%
      filter(p.value < pvalue_thresh)

    write_csv(counts_corr, peak_gene_outfile)
  }

main(
  metadata_infile = "${metadata_infile}",
  rna_counts_infile = "${rna_counts_infile}",
  atac_counts_infile = "${atac_counts_infile}",
  peak_gene_infile = "${peak_gene_infile}",
  peak_gene_outfile = "${peak_gene_outfile}",
  ncores = as.integer("${task.cpus}")
)
