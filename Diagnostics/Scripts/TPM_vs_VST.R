#Testing with TPM

library(GenomicFeatures)

library(rtracklayer)


# Path to your GFF3 file
gff_path <- "Montipora_capitata_HIv3.genes.gff3"

# Import GFF for ID mapping
gff <- import(gff_path)

# Create TxDb object from GFF
txdb <- makeTxDbFromGFF(gff_path, format = "gff3")


gff <- import("Montipora_capitata_HIv3.genes.gff3")

table(gff$type)

# Load GFF3 and build TxDb
txdb <- makeTxDbFromGFF("Montipora_capitata_HIv3.genes.gff3", format = "gff3")

# Step 2: Get exons grouped by transcript
exons_by_tx <- exonsBy(txdb, by = "tx")

reduced_exons <- GenomicRanges::reduce(exons_by_tx)

tx_lengths <- sum(width(reduced_exons))

colnames(mcols(transcripts))
#[1] "source" "type"   "score"  "phase"  "ID"     "Parent"
transcripts$transcript_id <- transcripts$ID

txdb_ids <- names(tx_lengths)

# Add proper transcript IDs (already present)
transcripts$transcript_id <- transcripts$ID

# Add txdb-style internal IDs (for matching tx_lengths)
transcripts$txdb_id <- as.character(seq_len(length(transcripts)))

# tx_lengths from earlier: named vector from exonsBy(txdb, by = "tx")
txdb_ids <- names(tx_lengths)

# Build transcript ID mapping
transcript_map <- data.frame(
  txdb_id = transcripts$txdb_id,
  transcript_id = transcripts$transcript_id
)

# Combine lengths with transcript IDs
tx_length_df <- data.frame(
  txdb_id = txdb_ids,
  length = as.numeric(tx_lengths)
)

# Merge the two
length_df <- merge(tx_length_df, transcript_map, by = "txdb_id")

# Final clean length table
length_df <- length_df[, c("transcript_id", "length")]

exons_tx <- exons_by_tx[["1"]]  # Or whatever the txdb ID is for that transcript
exons_tx


length_vector <- setNames(length_df$length, length_df$transcript_id)

matched_lengths <- length_vector[rownames(gcount_filt)]

# Optional: check for any missing IDs
sum(is.na(matched_lengths))  # Should be 0 ideally
#0 :)

counts_to_tpm <- function(counts, lengths) {
  rate <- counts / lengths
  tpm <- t(t(rate) / colSums(rate)) * 1e6
  return(tpm)
}

tpm_matrix <- counts_to_tpm(gcount_filt, matched_lengths)

log_tpm <- log2(tpm_matrix + 1)


# Calculate variance for each gene/transcript across samples
gene_vars <- apply(log_tpm, 1, var)

# Select top 500 (or 1000, or 2000)
top_genes <- head(order(gene_vars, decreasing = TRUE), 25741)

# Subset log_tpm
top_log_tpm <- log_tpm[top_genes, ]

# Transpose for PCA (samples as rows)
pca_res <- prcomp(t(top_log_tpm), scale. = TRUE)

pca_df <- as.data.frame(pca_res$x)
pca_df$sampleID <- rownames(pca_df)

# Join metadata
library(dplyr)
pca_df <- left_join(pca_df, treatmentinfo, by = "sampleID")

# Calculate % variance
percent_var <- round(100 * summary(pca_res)$importance[2, 1:2], 1)

library(ggplot2)
ggplot(pca_df, aes(PC1, PC2, color = timepoint, shape = timepoint)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  coord_fixed() +
  theme_classic() +
  ggtitle("PCA on log2(TPM + 1)")


tpm_pca <- prcomp(t(log_tpm), scale. = TRUE)

pca_df <- as.data.frame(tpm_pca$x)
pca_df$sampleID <- rownames(pca_df)

library(dplyr)
pca_df <- left_join(pca_df, treatmentinfo, by = "sampleID")

library(ggplot2)
ggplot(pca_df, aes(PC1, PC2, color = timepoint, shape = timepoint)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", round(summary(tpm_pca)$importance[2,1] * 100, 1), "% variance")) +
  ylab(paste0("PC2: ", round(summary(tpm_pca)$importance[2,2] * 100, 1), "% variance")) +
  theme_classic() +
  ggtitle("PCA on TPM (log2(TPM + 1))")

setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/P. acuta/Input")


# ─────────────────────────────────────────────────────────────────────
# Load necessary libraries
# ─────────────────────────────────────────────────────────────────────
library(GenomicFeatures)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(readr)
library(genefilter)

# ─────────────────────────────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────────────────────────────

# Load raw gene count matrix
gcount <- as.data.frame(read.csv("Pacu_gene_count_matrix_newGFF.csv", row.names = "gene_id"))

# Load sample metadata
treatmentinfo_pa <- read_csv("5-Pacu-SampleInfo.csv")

# ─────────────────────────────────────────────────────────────────────
# Filter low-expression transcripts
# ─────────────────────────────────────────────────────────────────────

# Filter: Keep genes with counts >10 in ≥33% of samples
filt <- filterfun(pOverA(0.33, 10))
gfilt <- genefilter(gcount, filt)

# Apply filter
gcount_filt_pa <- gcount[gfilt, ]

cat("Filtered from", nrow(gcount), "to", nrow(gcount_filt_pa), "transcripts.\n")

# ─────────────────────────────────────────────────────────────────────
# Load GFF and build transcript database
# ─────────────────────────────────────────────────────────────────────

# Load GFF3 file
gff_path <- "Pocillopora_acuta_HIv2.genes_fixed.gff3"
gff <- import(gff_path)

# Check available types
print(table(gff$type))  # Look for mRNA / transcript / gene

# Build TxDb
txdb <- makeTxDbFromGFF(gff_path, format = "gff3")

# ─────────────────────────────────────────────────────────────────────
# Get exon lengths per transcript
# ─────────────────────────────────────────────────────────────────────

exons_by_tx <- exonsBy(txdb, by = "tx")
reduced_exons <- GenomicRanges::reduce(exons_by_tx)
tx_lengths <- sum(width(reduced_exons))  # named vector

# ─────────────────────────────────────────────────────────────────────
# Extract transcript IDs from GFF
# ─────────────────────────────────────────────────────────────────────

# Use feature type that contains transcript IDs
transcripts <- gff[gff$type == "transcript"]

# Get transcript ID mapping
transcripts$transcript_id <- transcripts$ID
transcripts$txdb_id <- as.character(seq_len(length(transcripts)))

transcript_map <- data.frame(
  txdb_id = transcripts$txdb_id,
  transcript_id = transcripts$transcript_id
)

tx_length_df <- data.frame(
  txdb_id = names(tx_lengths),
  length = as.numeric(tx_lengths)
)

# Merge length and ID mapping
length_df <- merge(tx_length_df, transcript_map, by = "txdb_id")
length_df <- length_df[, c("transcript_id", "length")]

# ─────────────────────────────────────────────────────────────────────
# Match transcript lengths to count matrix
# ─────────────────────────────────────────────────────────────────────

length_vector_pa <- setNames(length_df$length, length_df$transcript_id)

# Check intersection
shared_ids <- intersect(rownames(gcount_filt_pa), transcripts$ID)
length(shared_ids)  # Should be a large number (ideally close to nrow(gcount_filt_pa))

# Subset both datasets to shared transcript IDs
gcount_matched <- gcount_filt_pa[shared_ids, ]
lengths_matched <- length_vector_pa[shared_ids]

# ─────────────────────────────────────────────────────────────────────
# TPM normalization
# ─────────────────────────────────────────────────────────────────────

counts_to_tpm <- function(counts, lengths) {
  rate <- counts / lengths
  tpm <- t(t(rate) / colSums(rate)) * 1e6
  return(tpm)
}

tpm_matrix <- counts_to_tpm(gcount_matched, lengths_matched)
log_tpm <- log2(tpm_matrix + 1)

# ─────────────────────────────────────────────────────────────────────
# PCA analysis
# ─────────────────────────────────────────────────────────────────────

# Calculate variance and select top variable genes (optional)
gene_vars <- apply(log_tpm, 1, var)
top_genes <- head(order(gene_vars, decreasing = TRUE), nrow(log_tpm))
top_log_tpm <- log_tpm[top_genes, ]

# PCA
pca_res <- prcomp(t(top_log_tpm), scale. = TRUE)
pca_df <- as.data.frame(pca_res$x)
pca_df$sampleID <- rownames(pca_df)

# Join metadata
pca_df <- left_join(pca_df, treatmentinfo_pa, by = "sampleID")

# % Variance
percent_var <- round(100 * summary(pca_res)$importance[2, 1:2], 1)

# Plot
p <- ggplot(pca_df, aes(PC1, PC2, color = timepoint, shape = timepoint)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  coord_fixed() +
  theme_classic() +
  ggtitle("PCA on Pocillopora acuta (log2(TPM + 1))")

print(p)


dim(log_tpm)
summary(pca_res)




# ─────────────────────────────────────────────────────────────
# Load Libraries
# ─────────────────────────────────────────────────────────────
library(DESeq2)
library(ggplot2)
library(dplyr)
library(readr)
library(genefilter)
library(patchwork)

# ─────────────────────────────────────────────────────────────
# Load Data (Pocillopora acuta)
# ─────────────────────────────────────────────────────────────
gcount_pa <- as.data.frame(read.csv("Pacu_gene_count_matrix_newGFF.csv", row.names = "gene_id"))
treatmentinfo_pa <- read_csv("5-Pacu-SampleInfo.csv")
treatmentinfo_pa$timepoint <- factor(treatmentinfo_pa$timepoint, levels = c("I", "II", "III"))

# ─────────────────────────────────────────────────────────────
# Filter Low-Count Genes
# ─────────────────────────────────────────────────────────────
filt_pa <- filterfun(pOverA(0.33, 10))
gfilt_pa <- genefilter(gcount_pa, filt_pa)
gcount_filt_pa <- gcount_pa[gfilt_pa, ]

cat("Genes before:", nrow(gcount_pa), "| Genes after filtering:", nrow(gcount_filt_pa), "\n")
#Genes before: 33730 | Genes after filtering: 20891 

# ─────────────────────────────────────────────────────────────
# VST PCA
# ─────────────────────────────────────────────────────────────
dds_pa <- DESeqDataSetFromMatrix(countData = gcount_filt_pa,
                                 colData = treatmentinfo_pa,
                                 design = ~ timepoint)

dds_pa <- estimateSizeFactors(dds_pa)
vst_data_pa <- vst(dds_pa, blind = FALSE)

vst_pca_df_pa <- plotPCA(vst_data_pa, intgroup = "timepoint", returnData = TRUE, ntop=20891)
percent_var_vst_pa <- round(100 * attr(vst_pca_df_pa, "percentVar"))

vst_pca_plot_pa <- ggplot(vst_pca_df_pa, aes(PC1, PC2, color = timepoint, shape = timepoint)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percent_var_vst_pa[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var_vst_pa[2], "% variance")) +
  ggtitle("PCA on VST-transformed counts") +
  coord_fixed() +
  theme_bw()


# ─────────────────────────────────────────────────────────────
# TPM PCA
# ─────────────────────────────────────────────────────────────
# Make sure your 'length_vector_pa' is available before this step
# It must be a named vector of transcript lengths (names = transcript IDs)

# Safety checks
shared_ids_pa <- intersect(rownames(gcount_filt_pa), names(length_vector_pa))
stopifnot(length(shared_ids_pa) > 0)
stopifnot(!any(is.na(length_vector_pa[shared_ids_pa])))

# Filter both to matching transcripts
gcount_matched_pa <- gcount_filt_pa[shared_ids_pa, ]
lengths_matched_pa <- length_vector_pa[shared_ids_pa]

# TPM calculation
counts_to_tpm <- function(counts, lengths) {
  rate <- counts / lengths
  tpm <- t(t(rate) / colSums(rate)) * 1e6
  return(tpm)
}

tpm_matrix_pa <- counts_to_tpm(gcount_matched_pa, lengths_matched_pa)
log_tpm_pa <- log2(tpm_matrix_pa + 1)

# PCA on TPM
pca_res_tpm_pa <- prcomp(t(log_tpm_pa), scale. = TRUE)
pca_df_tpm_pa <- as.data.frame(pca_res_tpm_pa$x)
pca_df_tpm_pa$sampleID <- rownames(pca_df_tpm_pa)
pca_df_tpm_pa <- left_join(pca_df_tpm_pa, treatmentinfo_pa, by = "sampleID")
percent_var_tpm_pa <- round(100 * summary(pca_res_tpm_pa)$importance[2, 1:2])

tpm_pca_plot_pa <- ggplot(pca_df_tpm_pa, aes(PC1, PC2, color = timepoint, shape = timepoint)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percent_var_tpm_pa[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var_tpm_pa[2], "% variance")) +
  ggtitle("PCA on log2(TPM + 1)") +
  coord_fixed() +
  theme_bw()


# ─────────────────────────────────────────────────────────────
# Combine and Display
# ─────────────────────────────────────────────────────────────
combined_pca_plot_pa <- vst_pca_plot_pa + tpm_pca_plot_pa + plot_layout(ncol = 2)
print(combined_pca_plot_pa)

# ─────────────────────────────────────────────────────────────
# Optional: Save
# ─────────────────────────────────────────────────────────────
ggsave("Pacu_VST_vs_TPM_PCA.png", combined_pca_plot_pa, width = 12, height = 6, dpi = 300)



