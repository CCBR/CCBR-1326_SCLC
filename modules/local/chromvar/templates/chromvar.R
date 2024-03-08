#!/usr/bin/env Rscript
# bioconductor packages
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
# tidyverse packages
library(dplyr)
library(readr)
library(stringr)
library(tibble)
# other cran packages
library(Matrix)

set.seed(2017)

# variables set via nextflow template
cpus <- as.integer("${task.cpus}") # 4
bed_filename <- "${consensus_bed}" # 'output/matrix_bed/raw_tmm_fpkm_batch_corrected_PDX_only.bed'
cluster_filename <- "${cluster_map}" # 'assets/cluster_membership.tsv'
output_tsv <- "${output_tsv}" # 'tmp.tsv'
output_rds <- "${output_rda}" # 'tmp.RData'

BiocParallel::register(BiocParallel::MulticoreParam(cpus, progressbar = TRUE))

cluster_dat <- read_tsv(cluster_filename) %>%
  mutate(filename = basename(bam))
bam_filenames <- cluster_dat %>% pull(filename)

peaks <- getPeaks(bed_filename, sort_peaks = TRUE)
seqinfo(peaks) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
peaks <- peaks %>%
  # resize can cause peaks to go over end of chromosomes
  GenomicRanges::resize(width = 500, fix = "center") %>%
  # trim down any peaks beyond end of chromosome
  GenomicRanges::trim()

# chromosomes with overflow
which(end(peaks) > seqlengths(BSgenome.Hsapiens.UCSC.hg19)[as.character(seqnames(peaks))])
GenomicRanges:::get_out_of_bound_index(peaks)

fragment_counts <- getCounts(bam_filenames,
  peaks,
  paired = TRUE,
  by_rg = FALSE,
  format = "bam",
  colData = DataFrame(ClusterName = cluster_dat %>% pull(ClusterName))
)

colSums(assay(fragment_counts))
fragment_counts

# verify there are no all zero rows
table(rowSums(assay(fragment_counts)) == 0)

rowData(fragment_counts)
fragment_counts <- addGCBias(fragment_counts,
  genome = BSgenome.Hsapiens.UCSC.hg19
)
rowData(fragment_counts)
hist(rowData(fragment_counts)[["bias"]])

# create background peaks with the same GC bias
bground <- getBackgroundPeaks(object = fragment_counts)

# colnames(rowData(example_counts))
# #find indices of samples to keep
# counts_filtered <- filterSamples(example_counts, min_depth = 1500,
#                                  min_in_peaks = 0.15, shiny = FALSE)
# counts_filtered <- filterPeaks(counts_filtered, non_overlapping = TRUE)
motifs <- getJasparMotifs()

# find motifs in ROI
motif_ix <- matchMotifs(motifs, fragment_counts,
  genome = BSgenome.Hsapiens.UCSC.hg19
)

# get deviations
dev <- computeDeviations(object = fragment_counts, annotations = motif_ix, background_peaks = bground)

# get variations
variability <- computeVariability(dev)

plotVariability(variability, use_plotly = FALSE)

vdf <- as.data.frame(variability)
devzdf <- as.data.frame(assays(dev)[["z"]])
rownames_to_column(vdf, var = "motif") -> vdf
vdf <- vdf[order(vdf[["p_value_adj"]], -vdf[["variability"]]), ]
rownames_to_column(devzdf, var = "motif") -> devzdf
df <- merge(vdf, devzdf, by = "motif")
write.table(df, file = output_tsv, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

topmotifs <- head(vdf["motif"], 100)
topdevzdf <- devzdf[(devzdf["motif"] %in% topmotifs), ]
rownames(topdevzdf) <- NULL
column_to_rownames(as.data.frame(topdevzdf), var = "motif") -> topdevzdf

save.image(output_rda)
