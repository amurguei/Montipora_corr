# --- Working dir (adjust if needed) ---
setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Montipora_correction/all_samples")

# if needed:
# install.packages(c("tidyverse","data.table","ggrepel"))
# if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
# BiocManager::install(c("genefilter","DESeq2"))

library(tidyverse)
library(data.table)
library(ggrepel)
library(genefilter)
library(DESeq2)

counts_file <- "counts_all.SRR.tsv"
meta_file   <- "metadata_aligned.to_counts_SRR.tsv"

# ---------- Load + enforce NUMERIC ----------
counts_df <- fread(counts_file)
gene_col  <- names(counts_df)[1]

# Build a numeric matrix explicitly
counts <- as.matrix(counts_df[, -1, with = FALSE])
rownames(counts) <- counts_df[[gene_col]]

# force numeric (character -> numeric); replace NA with 0
storage.mode(counts) <- "double"
counts[!is.finite(counts)] <- 0

# quick sanity check
stopifnot(is.matrix(counts), is.numeric(counts))
message(sprintf("Counts matrix: %d genes x %d samples", nrow(counts), ncol(counts)))

meta <- read.delim(meta_file, stringsAsFactors = FALSE)
stopifnot(all(meta$sample_id %in% colnames(counts)))
# align order
counts <- counts[, meta$sample_id, drop = FALSE]

# ---------- Plotting factors ----------
meta <- meta %>%
  mutate(
    BioProject = factor(batch),
    LifeStage  = factor(host_life_stage),
    hpf_label  = paste0(host_age_hpf, " hpf"),
    Species    = "OneSpecies"
  )

# ---------- FILTERING ----------
library(genefilter)

## starting point: counts_df -> gcount (numeric matrix)
# gcount <- as.matrix(counts_df[, -1])
# rownames(gcount) <- counts_df[[1]]
# storage.mode(gcount) <- "double"
# gcount[!is.finite(gcount)] <- 0

## pOverA filter: p = 3/48 = 0.0625, A = 10  (STRICT >10)
filt  <- filterfun(pOverA(0.0625, 10))
gfilt <- genefilter(gcount, filt)          # logical vector per gene

## Identify genes to keep by count filter
gkeep <- gcount[gfilt, , drop = FALSE]

## Identify gene list
gn.keep <- rownames(gkeep)

## Gene count data filtered in pOverA
gcount_filt <- as.data.frame(gcount[rownames(gcount) %in% gn.keep, , drop = FALSE])

## How many rows before and after?
cat("Before:", nrow(gcount), "\nAfter: ", nrow(gcount_filt), "\n")

## (optional) sanity — these should match
stopifnot(nrow(gkeep) == length(gn.keep),
          nrow(gkeep) == nrow(gcount_filt),
          sum(gfilt)  == nrow(gcount_filt))


library(tidyverse)
library(ggrepel)
library(DESeq2)

md <- meta_df

# add sample_id from rownames only if it's missing
if (!"sample_id" %in% names(md)) {
  md <- tibble::rownames_to_column(md, "sample_id")
}

# enforce uniqueness and type
# ---- deps (safe to rerun) ----
library(tidyverse)
library(ggrepel)
library(DESeq2)

# ---------- build a SAFE metadata table with exactly one 'sample_id' ----------
md <- meta_df
if (!"sample_id" %in% names(md)) {
  md <- tibble::rownames_to_column(md, "sample_id")
}
md <- md %>%
  mutate(sample_id = as.character(sample_id)) %>%
  distinct(sample_id, .keep_all = TRUE)

# align metadata to matrix columns (and check)
stopifnot(all(colnames(gmat) %in% md$sample_id))
md <- md[match(colnames(gmat), md$sample_id), ]
rownames(md) <- md$sample_id

# derive plotting fields (only if missing)
if (!"BioProject" %in% names(md)) md$BioProject <- factor(md$batch)
if (!"LifeStage"  %in% names(md)) md$LifeStage  <- factor(md$host_life_stage)
if (!"hpf_label"  %in% names(md)) md$hpf_label  <- paste0(md$host_age_hpf, " hpf")

# ---------- plotting helper (palettes from data actually present) ----------
make_pca_plot <- function(pca, df, title) {
  percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))[1:2]
  bp_levels <- levels(factor(df$BioProject))
  bp_pal <- setNames(
    c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")[seq_along(bp_levels)],
    bp_levels
  )
  ls_levels  <- levels(factor(df$LifeStage))
  shape_vals <- setNames(rep(21:25, length.out = length(ls_levels)), ls_levels)
  
  ggplot(df, aes(PC1, PC2)) +
    geom_point(aes(fill = BioProject, shape = LifeStage),
               size = 4.3, color = "black", stroke = 0.7) +
    ggrepel::geom_text_repel(aes(label = hpf_label),
                             size = 3, max.overlaps = 100,
                             box.padding = 0.35, point.padding = 0.3) +
    scale_fill_manual(values = bp_pal, name = "BioProject") +
    scale_shape_manual(values = shape_vals, name = "Life stage") +
    guides(
      fill  = guide_legend(override.aes = list(shape = 21, color = "black")),
      shape = guide_legend(override.aes = list(fill = "grey85", color = "black"))
    ) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(),
          legend.title = element_text(size = 11),
          legend.text  = element_text(size = 10),
          plot.title   = element_text(hjust = 0.5, face = "bold")) +
    ggtitle(title)
}

# =========================
# PCA 1: logCPM (CENTERED)
# =========================
libsize <- colSums(gmat)
cpm     <- sweep(gmat, 2, libsize/1e6, "/")
logcpm  <- log2(cpm + 1)

# <-- center = TRUE (key change)
pca1c <- prcomp(t(logcpm), center = TRUE, scale. = FALSE)
pca_df1c <- as_tibble(pca1c$x[, 1:2], rownames = "sample_id") %>%
  left_join(md, by = "sample_id")

p1c <- make_pca_plot(pca1c, pca_df1c, "PCA (logCPM; centered)")
plot(p1c)
ggsave("PCA_logCPM_centered.png", p1c, width = 7.2, height = 5.6, dpi = 600)

# =========================
# PCA 2: VST (design=~LifeStage, CENTERED)
# =========================
dds <- DESeqDataSetFromMatrix(countData = round(gmat),
                              colData   = DataFrame(md),
                              design    = ~ LifeStage)
dds <- dds[rowSums(counts(dds)) > 0, ]

vst_mat <- assay(vst(dds, blind = FALSE))

# <-- center = TRUE (key change)
pca2c <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)
pca_df2c <- as_tibble(pca2c$x[, 1:2], rownames = "sample_id") %>%
  left_join(md, by = "sample_id")

p2c <- make_pca_plot(pca2c, pca_df2c, "PCA (VST; LifeStage in design; centered)")

plot (p2c)
ggsave("PCA_VST_centered.png", p2c, width = 7.2, height = 5.6, dpi = 600)

# (optional) quick spread check
apply(pca1c$x[,1:2], 2, function(x) c(range=diff(range(x)), sd=sd(x)))
apply(pca2c$x[,1:2], 2, function(x) c(range=diff(range(x)), sd=sd(x)))

# ---- safety & prep ----
library(tidyverse)
library(ggrepel)


library(GenomicFeatures)

library(rtracklayer)
# Packages
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)

# Point to  GFF3 
gff_path <- "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/M. capitata/Inputs/Montipora_capitata_HIv3.genes.gff3"

# Build a TxDb
txdb <- makeTxDbFromGFF(gff_path, format = "gff3")

# EXON-UNION LENGTHS PER TRANSCRIPT  -------------------------
# use.names=TRUE tries to carry transcript names from the GFF (e.g., transcript_id/ID)
ebt <- exonsBy(txdb, by = "tx", use.names = TRUE)

# length per transcript = sum of reduced exon widths
tx_len <- vapply(ebt, function(gr) sum(width(reduce(gr))), integer(1))  # named integer vector

# (optional) if names are TXIDs and you need transcript names, map TXID -> TXNAME
if (!all(rownames(gcount_filt) %in% names(tx_len))) {
  k <- names(tx_len)                               # TXIDs
  map <- select(txdb, keys = k, keytype = "TXID", columns = c("TXID","TXNAME","GENEID"))
  txname_map <- setNames(map$TXNAME, map$TXID)
  # replace names when TXNAME exists, otherwise keep TXID
  new_names <- ifelse(!is.na(txname_map[k]), txname_map[k], k)
  names(tx_len) <- new_names
}

# 4) Match lengths to your filtered count matrix rows
gmat <- as.matrix(gcount_filt)                     # genes/transcripts x samples
storage.mode(gmat) <- "double"
# keep only entries with a positive finite length
lens <- tx_len[rownames(gmat)]
ok <- is.finite(lens) & lens > 0
gmat <- gmat[ok, , drop = FALSE]
lens <- lens[ok]

cat("Matched lengths:", sum(ok), "of", nrow(gcount_filt), "rows\n")
#Matched lengths: 24269 of 24269 rows
matched_lengths <- lens

# 5) TPM + log
counts_to_tpm <- function(counts, lengths) {
  rate <- counts / lengths
  t(t(rate) / colSums(rate)) * 1e6
}
tpm <- counts_to_tpm(gmat, lens)
log_tpm <- log2(tpm + 1)

# ---- Build safe metadata (one 'sample_id' col) ----
md <- meta_df

# ensure a single sample_id column
if (!"sample_id" %in% names(md)) {
  md <- tibble::rownames_to_column(md, "sample_id")
}
md$sample_id <- as.character(md$sample_id)
md <- dplyr::distinct(md, sample_id, .keep_all = TRUE)

# add plotting columns only if missing
if (!"BioProject" %in% names(md)) md$BioProject <- factor(md$batch)
if (!"LifeStage"  %in% names(md)) md$LifeStage  <- factor(md$host_life_stage)
if (!"hpf_label"  %in% names(md)) md$hpf_label  <- paste0(md$host_age_hpf, " hpf")

# align metadata to TPM columns
stopifnot(all(colnames(log_tpm) %in% md$sample_id))
md <- md[match(colnames(log_tpm), md$sample_id), ]
rownames(md) <- md$sample_id


# ---- PCA on TPM (CENTERED) ----
pca_tpm <- prcomp(t(log_tpm), center = TRUE, scale. = FALSE)
percentVar <- round(100 * (pca_tpm$sdev^2 / sum(pca_tpm$sdev^2)))[1:2]

# Join scores to metadata
pca_df <- tibble::as_tibble(pca_tpm$x[, 1:2], rownames = "sample_id") %>%
  dplyr::left_join(md, by = "sample_id")

# Palettes computed from data actually present
bp_levels <- levels(factor(pca_df$BioProject))
bp_pal <- setNames(
  c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")[seq_along(bp_levels)],
  bp_levels
)
ls_levels  <- levels(factor(pca_df$LifeStage))
shape_vals <- setNames(rep(21:25, length.out = length(ls_levels)), ls_levels)

# ---- Plot ----
library(ggplot2); library(ggrepel)

p_tpm <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(aes(fill = BioProject, shape = LifeStage),
             size = 4.3, color = "black", stroke = 0.7) +
  ggrepel::geom_text_repel(aes(label = hpf_label),
                           size = 3, max.overlaps = 100,
                           box.padding = 0.35, point.padding = 0.3) +
  scale_fill_manual(values = bp_pal, name = "BioProject") +
  scale_shape_manual(values = shape_vals, name = "Life stage") +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, color = "black")),
    shape = guide_legend(override.aes = list(fill = "grey85", color = "black"))
  ) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        legend.title = element_text(size = 11),
        legend.text  = element_text(size = 10),
        plot.title   = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("PCA (TPM; centered)")

print(p_tpm)
ggsave("PCA_TPM_centered.png", p_tpm, width = 7.2, height = 5.6, dpi = 600)

#ComBat-seq

# ====== ComBat-Seq (batch = BioProject, group = LifeStage) + VST + PCA ======
# deps
library(sva)
library(DESeq2)
library(tidyverse)
library(ggrepel)

# md already has: sample_id, batch/BioProject, host_life_stage (LifeStage), host_age_hpf, etc.
# Ensure fields exist and align to counts
stopifnot(all(colnames(gmat) %in% md$sample_id))
md <- md[match(colnames(gmat), md$sample_id), ]
rownames(md) <- md$sample_id
md$BioProject <- droplevels(factor(if ("BioProject" %in% names(md)) md$BioProject else md$batch))
md$LifeStage  <- droplevels(factor(md$LifeStage %||% md$host_life_stage))
md$LifeStage_model <- factor(make.names(as.character(md$LifeStage)))
md$hpf_num <- as.numeric(md$host_age_hpf)

# drop all-zero genes (faster, safer)
gmat <- as.matrix(gmat)
storage.mode(gmat) <- "double"
gmat <- gmat[rowSums(gmat) > 0, , drop = FALSE]

# ---- ComBat-Seq with continuous covariate (Hpf) ----
covar <- model.matrix(~ scale(hpf_num), data = md)  # preserves age trend without confounding
combat_counts <- ComBat_seq(
  counts    = round(gmat),
  batch     = md$BioProject,
  covar_mod = covar
)

# ---- VST (design-aware; keep LifeStage in design for downstream) ----
dds_c <- DESeqDataSetFromMatrix(countData = combat_counts,
                                colData   = DataFrame(md),
                                design    = ~ LifeStage_model)
dds_c <- dds_c[rowSums(counts(dds_c)) > 0, ]
vsd_c <- vst(dds_c, blind = FALSE)

# ---- PCA via plotPCA (return data to customize aesthetics) ----
pdat <- DESeq2::plotPCA(vsd_c,
                        intgroup = c("BioProject","LifeStage","hpf_label"),
                        returnData = TRUE,
                        ntop = nrow(vsd_c))   # use all filtered genes
percentVar <- round(100 * attr(pdat, "percentVar"))

bp_levels <- levels(factor(pdat$BioProject))
bp_pal <- setNames(c("#1b9e77","#d95f02","#7570b3","#e7298a")[seq_along(bp_levels)], bp_levels)
ls_levels  <- levels(factor(pdat$LifeStage))
shape_vals <- setNames(rep(21:25, length.out = length(ls_levels)), ls_levels)

plot_combatseq_vst <- ggplot(pdat, aes(PC1, PC2)) +
  geom_point(aes(fill = BioProject, shape = LifeStage),
             size = 4.3, color = "black", stroke = 0.7) +
  ggrepel::geom_text_repel(aes(label = hpf_label),
                           size = 3, max.overlaps = 100,
                           box.padding = 0.35, point.padding = 0.3) +
  scale_fill_manual(values = bp_pal, name = "BioProject") +
  scale_shape_manual(values = shape_vals, name = "Life stage") +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, color = "black")),
    shape = guide_legend(override.aes = list(fill = "grey85", color = "black"))
  ) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  coord_fixed() +
  theme_classic(base_size = 13) +
  ggtitle("ComBat-Seq (batch=BioProject, covar=Hpf) → VST PCA")

print(plot_combatseq_vst)
ggsave("PCA_ComBatSeq_VST_hpf.png", plot_combatseq_vst, width = 7.2, height = 5.6, dpi = 600)

#ComBat with covariates

library(dplyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(sva)

# ---------- Harmonize stages (keep metamorphosed types separate) ----------
md <- meta_df
if (!"sample_id" %in% names(md)) md <- tibble::rownames_to_column(md, "sample_id")
md <- md %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  mutate(
    BioProject = factor(if ("BioProject" %in% names(.)) BioProject else batch),
    LifeStage  = factor(if ("LifeStage"  %in% names(.)) LifeStage  else host_life_stage),
    hpf_label  = if (!"hpf_label" %in% names(.)) paste0(host_age_hpf, " hpf") else hpf_label
  )

# Merge only the ones you requested:
# - "Larvae" + "swimming larva" -> "Larva (swimming)"
# - "recruit" + "Attached Recruit" -> "Attached Recruit"
md$Stage6 <- recode(as.character(md$LifeStage),
                    "Larvae"              = "Larva (swimming)",
                    "swimming larva"      = "Larva (swimming)",
                    "recruit"             = "Attached Recruit",
                    .default = as.character(md$LifeStage)
)
# Order levels (6 distinct)
md$Stage6 <- factor(md$Stage6,
                    levels = c("Fertilized Egg","Embryo","Larva (swimming)",
                               "metamorphosed larva","Metamorphosed Polyp","Attached Recruit")
)

# Model-safe stage for DESeq2 design
md$Stage6_model <- factor(make.names(as.character(md$Stage6)))

# Align counts order
gmat <- as.matrix(gcount_filt)
storage.mode(gmat) <- "double"
stopifnot(all(colnames(gmat) %in% md$sample_id))
md <- md[match(colnames(gmat), md$sample_id), ]
rownames(md) <- md$sample_id

# ---------- ComBat-Seq with continuous Hpf covariate (avoids confounding) ----------
md$hpf_num <- as.numeric(md$host_age_hpf)
covar <- model.matrix(~ scale(hpf_num), data = md)

combat_counts <- ComBat_seq(
  counts    = round(gmat),
  batch     = droplevels(md$BioProject),
  covar_mod = covar
)

# ---------- VST (design-aware) ----------
dds_c <- DESeqDataSetFromMatrix(countData = combat_counts,
                                colData   = DataFrame(md),
                                design    = ~ Stage6_model)
dds_c <- dds_c[rowSums(counts(dds_c)) > 0, ]
vsd_c <- vst(dds_c, blind = FALSE)

# ---------- PCA via plotPCA ----------
pdat <- DESeq2::plotPCA(vsd_c,
                        intgroup = c("BioProject","Stage6","hpf_label"),
                        returnData = TRUE,
                        ntop = nrow(vsd_c))
percentVar <- round(100 * attr(pdat, "percentVar"))

# ---------- Aesthetics: color = BioProject (stroke), shape = Stage6 (6 distinct) ----------
# 6 distinct shapes for Stage6 (no repeats)
shape_map <- c(
  "Fertilized Egg"      = 21,  # filled circle
  "Embryo"              = 22,  # filled square
  "Larva (swimming)"    = 24,  # filled triangle up
  "metamorphosed larva" = 23,  # filled diamond
  "Metamorphosed Polyp" = 25,  # filled triangle down
  "Attached Recruit"    = 8    # outline star/asterisk-like (keeps it unique)
)
shapes_present <- shape_map[names(shape_map) %in% levels(factor(pdat$Stage6))]

bp_levels <- levels(factor(pdat$BioProject))
bp_pal <- setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")[seq_along(bp_levels)], bp_levels)

p_combatseq_vst_shapes <- ggplot(pdat, aes(PC1, PC2)) +
  geom_point(aes(color = BioProject, shape = Stage6),
             size = 4.6, stroke = 1.2) +
  ggrepel::geom_text_repel(aes(label = hpf_label),
                           size = 3, max.overlaps = 100,
                           box.padding = 0.35, point.padding = 0.3) +
  scale_color_manual(values = bp_pal, name = "BioProject") +
  scale_shape_manual(values = shapes_present, name = "Life stage") +
  guides(
    color = guide_legend(override.aes = list(shape = 21, size = 4.6, stroke = 1.0)),
    shape = guide_legend(override.aes = list(color = "grey20"))
  ) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  coord_fixed() +
  theme_classic(base_size = 13) +
  theme(
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  ) +
  ggtitle("Batch correction + VST")

print(p_combatseq_vst_shapes)
ggsave("PCA_ComBatSeq_VST_shapes.png", p_combatseq_vst_shapes, width = 7.6, height = 5.8, dpi = 600)

#Subsetting

# deps
library(DESeq2)
library(tidyverse)
library(ggrepel)

# ---- ensure harmonized stage labels exist (Stage6) ----
if (!"Stage6" %in% names(md)) {
  md$Stage6 <- dplyr::recode(as.character(md$LifeStage),
                             "Larvae"              = "Larva (swimming)",
                             "swimming larva"      = "Larva (swimming)",
                             "recruit"             = "Attached Recruit",
                             .default = as.character(md$LifeStage)
  )
}
md$Stage6 <- factor(md$Stage6,
                    levels = c("Fertilized Egg","Embryo","Larva (swimming)",
                               "metamorphosed larva","Metamorphosed Polyp","Attached Recruit")
)
if (!"hpf_label" %in% names(md)) md$hpf_label <- paste0(md$host_age_hpf, " hpf")
md$hpf_num <- as.numeric(md$host_age_hpf)

# ---- pick the three groups you asked for ----
keep_samp <- with(md,
                  (Stage6 == "Larva (swimming)"    & hpf_num == 163) |
                    (Stage6 == "metamorphosed larva") |
                    (Stage6 == "Attached Recruit"    & hpf_num == 231)
)
md_sub <- droplevels(md[keep_samp, ])

# align counts to subset
combat_sub <- combat_counts[, md_sub$sample_id, drop = FALSE]

# ---- VST on the subset (design-aware) ----
md_sub$Stage6_model <- factor(make.names(as.character(md_sub$Stage6)))
dds_sub <- DESeqDataSetFromMatrix(countData = combat_sub,
                                  colData   = DataFrame(md_sub),
                                  design    = ~ Stage6_model)
dds_sub <- dds_sub[rowSums(counts(dds_sub)) > 0, ]
vsd_sub <- vst(dds_sub, blind = FALSE)

# ---- PCA via plotPCA, then custom aesthetics ----
pdat <- DESeq2::plotPCA(vsd_sub,
                        intgroup = c("BioProject","Stage6","hpf_label"),
                        returnData = TRUE,
                        ntop = nrow(vsd_sub))
percentVar <- round(100 * attr(pdat, "percentVar"))

# shapes for the 3 stages (all distinct)
shape_map <- c("Larva (swimming)" = 24,  # triangle up
               "metamorphosed larva" = 23,  # diamond
               "Attached Recruit" = 22)     # square
shapes_present <- shape_map[names(shape_map) %in% levels(factor(pdat$Stage6))]

# colors by BioProject
bp_levels <- levels(factor(pdat$BioProject))
bp_pal <- setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")[seq_along(bp_levels)], bp_levels)

plot_three <- ggplot(pdat, aes(PC1, PC2)) +
  geom_point(aes(color = BioProject, shape = Stage6),
             size = 4.8, stroke = 1.1) +
  ggrepel::geom_text_repel(aes(label = hpf_label),
                           size = 3, max.overlaps = 100,
                           box.padding = 0.35, point.padding = 0.3) +
  scale_color_manual(values = bp_pal, name = "BioProject") +
  scale_shape_manual(values = shapes_present, name = "Life stage") +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  coord_fixed() +
  theme_classic(base_size = 13) +
  ggtitle("Subset + batch correction")

print(plot_three)
ggsave("PCA_ComBatSeq_VST_three_groups.png", plot_three, width = 7.4, height = 5.8, dpi = 600)


#With 280 hpf

library(DESeq2)
library(tidyverse)
library(ggrepel)

# Ensure harmonized stages exist (Stage6) + hpf numeric
if (!"Stage6" %in% names(md)) {
  md$Stage6 <- dplyr::recode(as.character(md$LifeStage),
                             "Larvae" = "Larva (swimming)", "swimming larva" = "Larva (swimming)",
                             "recruit" = "Attached Recruit", .default = as.character(md$LifeStage)
  )
}
md$hpf_num <- as.numeric(md$host_age_hpf)

# --- keep: 163 hpf Larva (swimming), ALL metamorphosed larva, 280 hpf Attached Recruit ---
keep_samp <- with(md,
                  (Stage6 == "Larva (swimming)"    & hpf_num == 163) |
                    (Stage6 == "metamorphosed larva") |
                    (Stage6 == "Attached Recruit"    & hpf_num == 280)
)
md_sub <- droplevels(md[keep_samp, ])
combat_sub <- combat_counts[, md_sub$sample_id, drop = FALSE]

# VST on subset
md_sub$Stage6_model <- factor(make.names(as.character(md_sub$Stage6)))
dds_sub <- DESeqDataSetFromMatrix(countData = combat_sub,
                                  colData   = DataFrame(md_sub),
                                  design    = ~ Stage6_model)
dds_sub <- dds_sub[rowSums(counts(dds_sub)) > 0, ]
vsd_sub <- vst(dds_sub, blind = FALSE)

# PCA + plot
pdat <- DESeq2::plotPCA(vsd_sub,
                        intgroup = c("BioProject","Stage6","hpf_label"),
                        returnData = TRUE,
                        ntop = nrow(vsd_sub))
percentVar <- round(100 * attr(pdat, "percentVar"))

shape_map <- c("Larva (swimming)" = 24, "metamorphosed larva" = 23, "Attached Recruit" = 22)
shapes_present <- shape_map[names(shape_map) %in% levels(factor(pdat$Stage6))]

bp_levels <- levels(factor(pdat$BioProject))
bp_pal <- setNames(c("#1b9e77","#d95f02","#7570b3","#e7298a")[seq_along(bp_levels)], bp_levels)

plot_three_280 <- ggplot(pdat, aes(PC1, PC2)) +
  geom_point(aes(color = BioProject, shape = Stage6), size = 4.8, stroke = 1.1) +
  ggrepel::geom_text_repel(aes(label = hpf_label), size = 3, max.overlaps = 100) +
  scale_color_manual(values = bp_pal, name = "BioProject") +
  scale_shape_manual(values = shapes_present, name = "Life stage") +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  coord_fixed() + theme_classic(base_size = 13) +
  ggtitle("160 hpf larvae, 230 recruit, 280 hpf Recruit)")

print(plot_three_280)
ggsave("Batch correction + VST + subset 280 hpf.png", plot_three_280, width = 7.4, height = 5.8, dpi = 600)
