#!/usr/bin/env Rscript
options(error = rlang::entrace)
library(dplyr)
library(furrr)
library(future)
library(glue)
library(purrr)
library(readr)
library(stringr)
library(tidyr)

main <-
  function(matched_ids_infile = "output_2/match_atac_rna_ids/matched_sample_ids.tsv",
           rna_counts_infile = "output_2/match_atac_rna_ids/rna_counts_norm_long.tsv",
           atac_counts_infile = "output_2/match_atac_rna_ids/atac_counts_norm_long.tsv",
           motifs_peaks_infile = "output_2/reformat_bed_intersect_motif/intersect.m334087_S35_mpbs.fixed.filt.intersect.raw_tmm_fpkm_batch_corrected_PDX_only.gencode.v19.annotation.TSS_padded.reformat.reformat_motif.tsv",
           outfile = "tmp.tsv",
           num_cores = 1) {
    motif_dat <- read_tsv(motifs_peaks_infile) %>%
      select(atac_sample_id, TF_gene_name, TSS_gene_name, peak_coord) %>%
      distinct()

    tf_lists <- motif_dat %>%
      # there can be multiple TFs per peak, so collapse them to a list
      group_by(atac_sample_id, peak_coord, TSS_gene_name) %>%
      summarize(TF_list = paste0(TF_gene_name, collapse = ","))

    motif_dat %>%
      # remove redundant TF-gene pairs
      select(-TF_gene_name) %>%
      distinct() %>%
      # join with TF lists
      left_join(tf_lists) %>%
      # join with matched atac-rna sample ids
      left_join(read_tsv(matched_ids_infile) %>%
        select(atac_sample_id, rna_sample_id, pdx_rank3) %>%
        rename(cluster = pdx_rank3)) %>%
      # join with atac counts
      left_join(
        read_tsv(atac_counts_infile) %>%
          rename(atac_sample_id = atac_id_sample) %>%
          select(-atac_id),
        by = c("atac_sample_id", "peak_coord")
      ) %>%
      # join with rna counts
      left_join(
        read_tsv(rna_counts_infile) %>%
          rename(rna_sample_id = rna_id_sample) %>%
          select(-rna_id),
        by = c("rna_sample_id", "TSS_gene_name" = "gene_name")
      ) %>%
      filter(!is.na(atac_count), !is.na(rna_count)) %>%
      write_tsv(outfile)
  }


args <- commandArgs(trailingOnly = TRUE)
# main()
main(
  matched_ids_infile = args[1],
  rna_counts_infile = args[2],
  atac_counts_infile = args[3],
  motifs_peaks_infile = args[4],
  outfile = args[5],
  num_cores = as.integer(args[6])
)
