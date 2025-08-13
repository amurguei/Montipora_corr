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

# ─────────────────────────────────────────────────────────────
# M. capitata (HIv3) — VST vs TPM PCA 
# ─────────────────────────────────────────────────────────────

  library(readr)
  library(dplyr)
  library(ggplot2)
  library(genefilter)
  library(DESeq2)
  library(patchwork)
  library(GenomicFeatures)
  library(AnnotationDbi)


options(stringsAsFactors = FALSE)

# ─────────────────────────────────────────────────────────────
# Paths & files
# ─────────────────────────────────────────────────────────────
wd  <- "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/M. capitata/Inputs"
out <- "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/TPM_PCA"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

counts_file     <- "Mcap_transcript_count_matrix.csv"
sampleinfo_file <- "5-Mcap-SampleInfo.csv"
gff_file        <- "Montipora_capitata_HIv3.genes.gff3"

stopifnot(file.exists(counts_file), file.exists(sampleinfo_file), file.exists(gff_file))

# ─────────────────────────────────────────────────────────────
# Load metadata & counts (masking-safe dplyr:: prefixes)
# ─────────────────────────────────────────────────────────────
tinfo      <- readr::read_csv(sampleinfo_file, show_col_types = FALSE)
gcount_raw <- readr::read_csv(counts_file, show_col_types = FALSE)

# Standardize ID column name
if ("Geneid" %in% names(gcount_raw)) gcount_raw <- dplyr::rename(gcount_raw, gene_id = Geneid)
if (!"gene_id" %in% names(gcount_raw) && "transcript_id" %in% names(gcount_raw)) {
  gcount_raw <- dplyr::rename(gcount_raw, gene_id = transcript_id)
}
stopifnot("gene_id" %in% names(gcount_raw))

# Keep only gene_id + sample columns present in metadata
tinfo$sampleID <- as.character(tinfo$sampleID)
gcount <- gcount_raw |>
  dplyr::select(dplyr::any_of(c("gene_id", tinfo$sampleID))) |>
  as.data.frame()

missing_cols <- setdiff(tinfo$sampleID, colnames(gcount))
if (length(missing_cols)) stop("Missing sample columns in counts: ", paste(missing_cols, collapse = ", "))

rownames(gcount) <- gcount$gene_id
gcount <- gcount |>
  dplyr::mutate(dplyr::across(-gene_id, ~as.integer(round(.)))) |>
  dplyr::select(-gene_id)

# Order metadata to match count columns; ensure timepoint factor
tinfo <- dplyr::left_join(
  tibble::tibble(sampleID = colnames(gcount)),
  dplyr::select(tinfo, sampleID, timepoint),
  by = "sampleID"
)
stopifnot(identical(tinfo$sampleID, colnames(gcount)))
stopifnot("timepoint" %in% names(tinfo))
tinfo$timepoint <- factor(tinfo$timepoint, levels = c("I","II","III"))

# ─────────────────────────────────────────────────────────────
# Build transcript lengths from GFF
# ─────────────────────────────────────────────────────────────
txdb <- GenomicFeatures::makeTxDbFromGFF(gff_file, format = "gff3")
ex_by_tx <- exonsBy(txdb, by = "tx")
tx_lengths <- sum(width(reduce(ex_by_tx)))  # reduced exons per transcript

# Map TXID -> TXNAME to get public transcript IDs
tx_map <- AnnotationDbi::select(txdb, keys = names(tx_lengths),
                                keytype = "TXID", columns = "TXNAME")
txname_by_txid <- setNames(tx_map$TXNAME, tx_map$TXID)
names(tx_lengths) <- unname(txname_by_txid[names(tx_lengths)])
tx_lengths <- tx_lengths[!is.na(names(tx_lengths))]

length_vector <- setNames(as.numeric(tx_lengths), as.character(names(tx_lengths)))

# Sanity: overlap between counts and transcript IDs from GFF
prop_overlap <- mean(rownames(gcount) %in% names(length_vector))
message(sprintf("Transcript ID overlap with GFF: %.1f%%", 100*prop_overlap))
if (prop_overlap < 0.5) stop("Low overlap: are your counts gene-level? If so, switch to gene lengths.")

# ─────────────────────────────────────────────────────────────
# Filter low-count features (P=0.33 over A=10)
# ─────────────────────────────────────────────────────────────
filt <- filterfun(pOverA(0.33, 10))
gfilt <- genefilter(gcount, filt)
gcount_filt <- gcount[gfilt, , drop = FALSE]
cat("Features before:", nrow(gcount), "| after filtering:", nrow(gcount_filt), "\n")

# ─────────────────────────────────────────────────────────────
# VST PCA
# ─────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(gcount_filt)),
                              colData   = as.data.frame(tinfo),
                              design    = ~ timepoint)
dds <- estimateSizeFactors(dds)
vst_data <- vst(dds, blind = FALSE)

vst_df <- plotPCA(vst_data, intgroup = "timepoint", returnData = TRUE, ntop = nrow(gcount_filt))
pv_vst <- round(100 * attr(vst_df, "percentVar"))

p_vst <- ggplot(vst_df, aes(PC1, PC2, color = timepoint, shape = timepoint)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", pv_vst[1], "% variance")) +
  ylab(paste0("PC2: ", pv_vst[2], "% variance")) +
  ggtitle("PCA on VST-transformed counts") +
  coord_fixed() +
  theme_bw()

# ─────────────────────────────────────────────────────────────
# TPM PCA
# ─────────────────────────────────────────────────────────────
shared_ids <- intersect(rownames(gcount_filt), names(length_vector))
stopifnot(length(shared_ids) > 0, !any(is.na(length_vector[shared_ids])))

counts_to_tpm <- function(counts, lengths) {
  rate <- sweep(as.matrix(counts), 1, lengths, "/")
  t(t(rate) / colSums(rate)) * 1e6
}

tpm_mat  <- counts_to_tpm(gcount_filt[shared_ids, , drop = FALSE], length_vector[shared_ids])
log_tpm  <- log2(tpm_mat + 1)

pca_tpm  <- prcomp(t(log_tpm), scale. = TRUE)
df_tpm   <- as.data.frame(pca_tpm$x)
df_tpm$sampleID <- rownames(df_tpm)
df_tpm   <- dplyr::left_join(df_tpm, tinfo, by = "sampleID")
pv_tpm   <- round(100 * summary(pca_tpm)$importance[2, 1:2])

p_tpm <- ggplot(df_tpm, aes(PC1, PC2, color = timepoint, shape = timepoint)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", pv_tpm[1], "% variance")) +
  ylab(paste0("PC2: ", pv_tpm[2], "% variance")) +
  ggtitle("PCA on log2(TPM + 1)") +
  coord_fixed() +
  theme_bw()

# ─────────────────────────────────────────────────────────────
# Combine & save
# ─────────────────────────────────────────────────────────────
combined <- p_vst + p_tpm + patchwork::plot_layout(ncol = 2)
print(combined)

outfile <- file.path(out, "Mcap_VST_vs_TPM_PCA.png")
ggsave(outfile, combined, width = 12, height = 6, dpi = 300)
message("Saved: ", outfile)


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

# ─────────────────────────────────────────────────────────────
# Load Libraries
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(genefilter)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

# ─────────────────────────────────────────────────────────────
# Stylophora pistillata
# ─────────────────────────────────────────────────────────────
wd  <- "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/S. pistillata/Transcriptomics/Inputs"
out <- "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/TPM_PCA"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

species_label <- "S. pistillata"
species_tag   <- "Spis"

# ─────────────────────────────────────────────────────────────
# Filenames
# ─────────────────────────────────────────────────────────────
counts_file     <- "4-Spis-GeneCountMatrix.csv"   # count matrix
sampleinfo_file <- "5-Spis-SampleInfo.csv"        # metadata
lengths_file    <- "averaged_gene_lengths.csv"    # gene lengths

rowname_col <- "gene_id"

# ─────────────────────────────────────────────────────────────
# Load Data
# ─────────────────────────────────────────────────────────────
if (!file.exists(counts_file)) stop("Counts file not found: ", counts_file)
if (!file.exists(sampleinfo_file)) stop("Sample info file not found: ", sampleinfo_file)
if (!file.exists(lengths_file)) stop("Lengths file not found: ", lengths_file)

gcount <- as.data.frame(read.csv(counts_file, row.names = rowname_col, check.names = FALSE))
tinfo  <- read_csv(sampleinfo_file, show_col_types = FALSE)

# Load lengths
len_df <- read_csv(lengths_file, show_col_types = FALSE)
names(len_df) <- tolower(names(len_df))
length_vector <- as.numeric(len_df$average_length)
names(length_vector) <- as.character(len_df$gene_id)

# ─────────────────────────────────────────────────────────────
# Match samples
# ─────────────────────────────────────────────────────────────
common_samples <- intersect(colnames(gcount), tinfo$sampleID)
if (length(common_samples) == 0) stop("No overlapping samples between counts and sample info.")
gcount <- gcount[, common_samples, drop = FALSE]
tinfo  <- tinfo %>% filter(sampleID %in% common_samples)
tinfo  <- tinfo[match(colnames(gcount), tinfo$sampleID), , drop = FALSE]

# Ensure timepoint factor
tinfo$timepoint <- factor(tinfo$timepoint, levels = c("I", "II", "III"))

# ─────────────────────────────────────────────────────────────
# Filter Low-Count Genes
# ─────────────────────────────────────────────────────────────
filt  <- filterfun(pOverA(0.33, 10))
gfilt <- genefilter(gcount, filt)
gcount_filt <- gcount[gfilt, , drop = FALSE]

cat("Genes before:", nrow(gcount), "| Genes after filtering:", nrow(gcount_filt), "\n")

# ─────────────────────────────────────────────────────────────
# VST PCA
# ─────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(gcount_filt)),
                              colData   = as.data.frame(tinfo),
                              design    = ~ timepoint)
dds <- estimateSizeFactors(dds)
vst_data <- vst(dds, blind = FALSE)

vst_pca_df <- plotPCA(vst_data, intgroup = "timepoint", returnData = TRUE, ntop = nrow(gcount_filt))
percent_var_vst <- round(100 * attr(vst_pca_df, "percentVar"))

vst_pca_plot <- ggplot(vst_pca_df, aes(PC1, PC2, color = timepoint, shape = timepoint)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percent_var_vst[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var_vst[2], "% variance")) +
  ggtitle(paste0(species_label, " — PCA on VST-transformed counts")) +
  coord_fixed() +
  theme_bw()

# ─────────────────────────────────────────────────────────────
# TPM PCA
# ─────────────────────────────────────────────────────────────
shared_ids <- intersect(rownames(gcount_filt), names(length_vector))
if (length(shared_ids) == 0) stop("No overlap between gene IDs in counts and length vector.")
if (any(is.na(length_vector[shared_ids]))) stop("NA values in matched length vector.")

gcount_matched  <- gcount_filt[shared_ids, , drop = FALSE]
lengths_matched <- length_vector[shared_ids]

counts_to_tpm <- function(counts, lengths) {
  rate <- sweep(as.matrix(counts), 1, lengths, "/")
  t(t(rate) / colSums(rate)) * 1e6
}

tpm_matrix <- counts_to_tpm(gcount_matched, lengths_matched)
log_tpm    <- log2(tpm_matrix + 1)

pca_res_tpm <- prcomp(t(log_tpm), scale. = TRUE)
pca_df_tpm  <- as.data.frame(pca_res_tpm$x)
pca_df_tpm$sampleID <- rownames(pca_df_tpm)
pca_df_tpm  <- left_join(pca_df_tpm, tinfo, by = "sampleID")
percent_var_tpm <- round(100 * summary(pca_res_tpm)$importance[2, 1:2])

tpm_pca_plot <- ggplot(pca_df_tpm, aes(PC1, PC2, color = timepoint, shape = timepoint)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percent_var_tpm[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var_tpm[2], "% variance")) +
  ggtitle(paste0(species_label, " — PCA on log2(TPM + 1)")) +
  coord_fixed() +
  theme_bw()

print(tpm_pca_plot)

# ─────────────────────────────────────────────────────────────
# Combine & Save
# ─────────────────────────────────────────────────────────────
combined_pca_plot <- vst_pca_plot + tpm_pca_plot + plot_layout(ncol = 2)
print(combined_pca_plot)

outfile <- file.path(out, paste0(species_tag, "_VST_vs_TPM_PCA.png"))
ggsave(outfile, combined_pca_plot, width = 12, height = 6, dpi = 300)
message("Saved: ", outfile)


# ─────────────────────────────────────────────────────────────
# Acropora tenuis (Cooke / ReefGenomics) — PCA pipeline
# ─────────────────────────────────────────────────────────────
  library(tidyverse)
  library(genefilter)
  library(DESeq2)
  library(ggrepel)
  library(patchwork)
  library(GenomicFeatures)


options(stringsAsFactors = FALSE)

# ─────────────────────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────────────────────
wd  <- "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/A_tenuis"
out <- "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/TPM_PCA"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

species_label <- "Acropora tenuis (Cooke GBR)"
species_tag   <- "Aten"

# ─────────────────────────────────────────────────────────────
# Input files (edit these two if names differ)
# ─────────────────────────────────────────────────────────────
counts_file     <- "gene_count_matrix_Acropora.csv"  # e.g., your counts matrix
sampleinfo_file <- "Sample_info_A_tenuis.csv"        # e.g., has columns: sampleID, timepoint
gff_file        <- "aten_0.11.maker_post_001.genes.gff"

stopifnot(file.exists(counts_file), file.exists(sampleinfo_file), file.exists(gff_file))

# ─────────────────────────────────────────────────────────────
# Build/Load gene lengths from GFF (once)
# ─────────────────────────────────────────────────────────────
length_csv <- "averaged_gene_lengths_Atenuis.csv"
if (!file.exists(length_csv)) {
  message("Building gene-length table from GFF (first run only)…")
  txdb <- makeTxDbFromGFF(gff_file, format = "gff3")
  ex_by_gene <- exonsBy(txdb, by = "gene")
  # Reduce to avoid double-counting overlapping exons across isoforms
  gene_lengths <- sum(width(reduce(ex_by_gene)))
  length_df <- tibble(
    gene_id = names(gene_lengths),
    average_length = as.numeric(gene_lengths)
  )
  write_csv(length_df, length_csv)
  message("Saved gene lengths → ", normalizePath(length_csv))
} else {
  message("Using cached gene lengths: ", normalizePath(length_csv))
}

# ─────────────────────────────────────────────────────────────
# Load counts + sample info
# ─────────────────────────────────────────────────────────────
treatmentinfo <- read_csv(sampleinfo_file, show_col_types = FALSE)
gcount_raw    <- read_csv(counts_file, show_col_types = FALSE)

# Standardize counts format
if ("Geneid" %in% names(gcount_raw)) names(gcount_raw)[names(gcount_raw) == "Geneid"] <- "gene_id"
stopifnot("gene_id" %in% names(gcount_raw))

# keep only gene_id + sample columns (from metadata)
if ("Geneid" %in% names(gcount_raw)) gcount_raw <- dplyr::rename(gcount_raw, gene_id = Geneid)

gcount <- gcount_raw %>%
  dplyr::select(dplyr::any_of(c("gene_id", treatmentinfo$sampleID))) %>%
  as.data.frame()

rownames(gcount) <- gcount$gene_id
gcount <- gcount %>%
  dplyr::mutate(dplyr::across(-gene_id, as.integer)) %>%
  dplyr::select(-gene_id)

# make sure IDs are characters
treatmentinfo$sampleID <- as.character(treatmentinfo$sampleID)

treatmentinfo <- treatmentinfo %>%
  dplyr::filter(sampleID %in% colnames(gcount)) %>%
  dplyr::slice(match(colnames(gcount), sampleID))

# sanity
missing <- setdiff(colnames(gcount), treatmentinfo$sampleID)
if (length(missing)) stop("Missing in sample info: ", paste(missing, collapse=", "))

# Quick sanity: GFF ↔ counts gene ID overlap
len_df <- read_csv(length_csv, show_col_types = FALSE)
length_vector <- setNames(as.numeric(len_df$average_length), as.character(len_df$gene_id))
prop_overlap <- mean(rownames(gcount) %in% names(length_vector))
message(sprintf("Gene ID overlap with GFF: %.1f%%", 100*prop_overlap))

# ─────────────────────────────────────────────────────────────
# Filter low-count genes (P=0.33 over A=10)
# ─────────────────────────────────────────────────────────────
filt         <- filterfun(pOverA(0.33, 10))
gfilt        <- genefilter(gcount, filt)
gcount_filt  <- gcount[gfilt, , drop = FALSE]
cat("Genes before:", nrow(gcount), "| after filtering:", nrow(gcount_filt), "\n")

# ─────────────────────────────────────────────────────────────
# VST PCA (all samples, for inspection)
# ─────────────────────────────────────────────────────────────
gdds <- DESeqDataSetFromMatrix(countData = round(as.matrix(gcount_filt)),
                               colData   = as.data.frame(treatmentinfo),
                               design    = ~ timepoint)
gdds <- estimateSizeFactors(gdds)
gvst <- vst(gdds, blind = FALSE)

gPCAdata <- plotPCA(gvst, intgroup = "timepoint", returnData = TRUE, ntop = nrow(gcount_filt))
percentVar <- round(100 * attr(gPCAdata, "percentVar"))

p_all <- ggplot(gPCAdata, aes(PC1, PC2)) +
  geom_point(aes(shape = timepoint, colour = timepoint), size = 6) +
  ggrepel::geom_text_repel(aes(label = rownames(gPCAdata)), size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + theme_classic() + ggtitle(paste0(species_label, " — VST PCA (all samples)"))

print (p_all)
# ─────────────────────────────────────────────────────────────
# Remove unwanted samples 
# ─────────────────────────────────────────────────────────────
remove_samples <- c("DRR318292", "DRR318296", "DRR318290")  

treatmentinfo_f <- treatmentinfo[!treatmentinfo$sampleID %in% remove_samples, , drop = FALSE]
gvst_f          <- gvst[, !colnames(gvst) %in% remove_samples]

# Align and sanity-check
datExpr <- as.data.frame(t(assay(gvst_f)))
treatmentinfo_f <- treatmentinfo_f[match(rownames(datExpr), treatmentinfo_f$sampleID), , drop = FALSE]
stopifnot(all(rownames(datExpr) == treatmentinfo_f$sampleID))

gPCAdata_f <- plotPCA(gvst_f, intgroup = "timepoint", returnData = TRUE, ntop = ncol(datExpr))
percentVar_f <- round(100 * attr(gPCAdata_f, "percentVar"))

p_post <- ggplot(gPCAdata_f, aes(PC1, PC2)) +
  geom_point(aes(shape = timepoint, colour = timepoint), size = 6) +
  xlab(paste0("PC1: ", percentVar_f[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_f[2], "% variance")) +
  coord_fixed() + theme_classic() + ggtitle("VST PCA")

# ─────────────────────────────────────────────────────────────
# TPM PCA (after sample removal)
# ─────────────────────────────────────────────────────────────
counts_to_tpm <- function(counts, lengths) {
  rate <- sweep(as.matrix(counts), 1, lengths, "/")
  t(t(rate) / colSums(rate)) * 1e6
}

# Use raw counts (filtered genes), restricted to kept samples
gcount_post <- gcount_filt[, colnames(gvst_f), drop = FALSE]
shared_ids  <- intersect(rownames(gcount_post), names(length_vector))
stopifnot(length(shared_ids) > 0, !any(is.na(length_vector[shared_ids])))

tpm      <- counts_to_tpm(gcount_post[shared_ids, , drop = FALSE], length_vector[shared_ids])
log_tpm  <- log2(tpm + 1)
pca_tpm  <- prcomp(t(log_tpm), scale. = TRUE)
df_tpm   <- as.data.frame(pca_tpm$x) %>% rownames_to_column("sampleID") %>% left_join(as.data.frame(treatmentinfo_f), by = "sampleID")
pv_tpm   <- round(100 * summary(pca_tpm)$importance[2, 1:2])

p_tpm <- ggplot(df_tpm, aes(PC1, PC2)) +
  geom_point(aes(shape = timepoint, colour = timepoint), size = 6) +
  xlab(paste0("PC1: ", pv_tpm[1], "% variance")) +
  ylab(paste0("PC2: ", pv_tpm[2], "% variance")) +
  coord_fixed() + theme_classic() + ggtitle("TPM PCA")

plot(p_tpm)
# ─────────────────────────────────────────────────────────────
# Save combined figure
# ─────────────────────────────────────────────────────────────
combo <- (p_post + p_tpm)
print(combo)

outfile <- file.path(out, paste0(species_tag, "Acropora_tenuis_VST_TPM.png"))
ggsave(outfile, combo, width = 14, height = 12, dpi = 300)
message("Saved figure: ", outfile)



