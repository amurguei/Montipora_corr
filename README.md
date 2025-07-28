# Montipora_correction
Attempts to correct artifacts in Montipora capitata developmental timeseries (3 life stages)

The *Montipora capitata* developmental time series appears to contain a technical artifact inflating separation between the spat samples and other life stages. This artifact may influence cross-species comparative analyses.

### Remapping and Quality Control

The original mappings were done on Andromeda (now offline), and the MultiQC reports were only on raw and cleaned reads, to get more informative statistics, I:
- I **re-mapped all samples using Hive**.
- I included **Hisat2 logs**, **featureCounts summaries** (only for diagnostic purposes, I had previously used the prepDE script), and computed **Transcript Integrity Numbers (TINs)**.
- These QC metrics showed that **spat samples had lower alignment rates and TINs**, suggesting RNA degradation.

###  Correction Strategies

To reduce the technical separation, I tested three groups of corrections and visualized results via PCA:

1. **Incorporating Covariates into VST**
   - I included covariates in the DESeq2 model used for `vst()`:
     - Alignment rate
     - TINs
     - Surrogate variables from `svaseq()`

2. **Residualization of covariates**
   - I regressed out the same covariates out of the expression matrix, this created residualized matrices.

3. **ComBat-seq / Hypothetical Batch Adjustment**
   - I created a `spat_status` variable assuming possible batch effects (e.g., different sequencing run) that could help explain the observed pattern, assuming 'spat' belonged to a different batch as 'larvae' and 'metamorphosed'. 

### Results from PCA

- **VST with covariates**: No dramatic change in variance explained by PC1/PC2, but **spread along PC1 (range) shrank**, especially with SVA.
- **Residualized matrices**: These reduced variance but often distorted structure.
- **ComBat-seq**: Also introduced irregularities in global structure.

### WGCNA Network Comparison

I tested whether including SVs as covariates led to changes in the network. 

- I compared WGCNA networks built on:
  1. The **uncorrected VST matrix**
  2. The **VST matrix with SVA as a covariate**
- Results:
  - Core network structure remained **largely preserved**.
  - SVA correction slightly **redistributed module membership** (e.g., fewer genes in turquoise, more in smaller modules).

In this repository I placed relevant files for the mapping process (MultiQC reports, bioinformatics code), the WGCNA Montipora code and gene count matrix, and the attempts at corrections code with files. # Montipora_correction

