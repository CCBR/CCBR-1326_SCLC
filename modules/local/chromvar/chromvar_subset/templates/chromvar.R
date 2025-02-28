#!/usr/bin/env Rscript
# bioconductor packages
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
# tidyverse packages
library(dplyr)
library(ggplot2)
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
output_rda <- "${output_rda}" # 'tmp.RData'
output_png <- "${output_png}"
cluster_id <- "${cluster_id}"

BiocParallel::register(BiocParallel::MulticoreParam(cpus, progressbar = TRUE))

cluster_dat <- read_tsv(cluster_filename) %>%
  mutate(filename = basename(bam)) %>%
  filter(clusterName == cluster_id)
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
  colData = DataFrame(clusterName = cluster_dat %>% pull(clusterName))
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
# adapted from https://github.com/GreenleafLab/chromVAR/blob/0f27fcc8463d537d770f164a55f949701beb6add/R/motifs.R
getJasparMotifs_ccbr <- function(species = "Homo sapiens",
                                 collection = "CORE",
                                 jaspar_db = JASPAR2016::JASPAR2016,
                                 ...) {
  opts <- list()
  opts["species"] <- species
  opts["collection"] <- collection
  opts <- c(opts, list(...))
  out <- TFBSTools::getMatrixSet(jaspar_db, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out)))) {
    names(out) <- paste(names(out), TFBSTools::name(out), sep = "_")
  }
  return(out)
}
motifs <- getJasparMotifs_ccbr(jaspar_db = JASPAR2020::JASPAR2020)


# find motifs in ROI
motif_ix <- matchMotifs(motifs, fragment_counts,
  genome = BSgenome.Hsapiens.UCSC.hg19
)

# get deviations
dev <- computeDeviations(object = fragment_counts, annotations = motif_ix, background_peaks = bground)

# get variations
variability <- computeVariability(dev)

var_plot <- plotVariability(variability, use_plotly = FALSE)
ggsave(filename = output_png, plot = var_plot)

n_top <- 20
top_tfs <- variability %>%
  slice_max(order_by = variability, n = n_top) %>%
  arrange(desc(variability)) %>%
  mutate(rank = row_number()) %>%
  arrange(variability)
var_plot_top <- top_tfs %>%
  mutate(tf = factor(name, levels = top_tfs %>% pull(name))) %>%
  ggplot(aes(
    y = tf,
    x = variability,
    xmin = bootstrap_lower_bound,
    xmax = bootstrap_upper_bound
  )) +
  geom_pointrange() +
  labs(
    y = "",
    title = glue::glue("Top {n_top} TFs by variability in accessibility")
  ) +
  theme_bw()

tsne_results <- deviationsTsne(dev, threshold = 1.5, perplexity = 10)
tsne_plots <- plotDeviationsTsne(dev, tsne_results,
  annotation_name = "POU5F1B",
  sample_column = "clusterName",
  shiny = FALSE
)

vdf <- as.data.frame(variability)
devzdf <- as.data.frame(assays(dev)[["z"]])
rownames_to_column(vdf, var = "motif") -> vdf
vdf <- vdf[order(vdf[["p_value_adj"]], -vdf[["variability"]]), ]
rownames_to_column(devzdf, var = "motif") -> devzdf
df <- merge(vdf, devzdf, by = "motif")
write.table(df, file = output_tsv, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

topmotifs <- head(vdf[["motif"]], 100)
topdevzdf <- devzdf[(devzdf[["motif"]] %in% topmotifs), ]
rownames(topdevzdf) <- NULL
column_to_rownames(as.data.frame(topdevzdf), var = "motif") -> topdevzdf

save.image(output_rda)
