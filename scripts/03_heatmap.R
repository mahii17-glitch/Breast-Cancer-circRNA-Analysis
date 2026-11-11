# GSE192849: Reproducible circRNA heatmap
# Run from the Breast_Cancer_RNAseq_Project directory

if (!requireNamespace('pheatmap', quietly = TRUE)) stop('Install pheatmap first')

counts <- read.csv('data/cleaned_counts_matrix.csv', row.names = 1, check.names = FALSE)
sig <- read.csv('results/significant_circRNAs.csv', row.names = 1)
stopifnot(all(rownames(sig) %in% rownames(counts)))
selected <- as.matrix(counts[rownames(sig), , drop = FALSE])

size_factors <- c(Tumor1 = 1.982503, Normal1 = 0.510856,
                  Tumor2 = 0.630649, Normal2 = 1.822129,
                  Tumor3 = 1.194504, Normal3 = 0.994507)
stopifnot(all(colnames(selected) %in% names(size_factors)))

normalized <- sweep(selected, 2, size_factors[colnames(selected)], '/')
log_counts <- log2(normalized + 1)
stopifnot(all(apply(log_counts, 1, sd) > 0))
heatmap_z <- t(scale(t(log_counts)))

annotation <- data.frame(
  Patient = factor(rep(paste0('Patient', 1:3), each = 2)),
  Condition = factor(rep(c('Tumor', 'Normal'), 3))
)
rownames(annotation) <- colnames(heatmap_z)

dir.create('plots', showWarnings = FALSE)
pheatmap::pheatmap(
  heatmap_z, scale = 'none',
  cluster_rows = TRUE, cluster_cols = FALSE,
  annotation_col = annotation,
  main = 'Differentially Expressed circRNAs in Breast Cancer',
  fontsize_row = 8, fontsize_col = 10, angle_col = 45,
  border_color = NA,
  filename = 'plots/circRNA_heatmap_reproducible.png',
  width = 11, height = 8
)
