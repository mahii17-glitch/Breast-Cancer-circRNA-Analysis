# GSE192849: Paired circRNA differential expression
# Run 01_data_loading.R first

library(DESeq2)

# Filter low-count circRNAs
counts_filtered <- counts_int[
  rowSums(counts_int) >= 10, , drop = FALSE
]

# Sample metadata
sample_info <- data.frame(
  row.names = colnames(counts_filtered),
  condition = factor(
    c('Tumor','Normal','Tumor','Normal','Tumor','Normal'),
    levels = c('Normal','Tumor')
  ),
  patient = factor(rep(paste0('Patient', 1:3), each = 2))
)

# Create paired DESeq2 model
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = sample_info,
  design = ~ patient + condition
)

# Run differential expression analysis
dds <- DESeq(dds)

res <- as.data.frame(results(
  dds,
  contrast = c('condition','Tumor','Normal')
))

# Select significant circRNAs
sig_results <- res[
  !is.na(res$padj) & res$padj < 0.05,
  , drop = FALSE
]

sig_results <- sig_results[
  order(sig_results$padj), , drop = FALSE
]

cat('Filtered circRNAs:', nrow(counts_filtered), '\n')
cat('Significant circRNAs:', nrow(sig_results), '\n')

# Save differential expression results
dir.create('results', showWarnings = FALSE)
write.csv(res, 'results/all_differential_expression_results.csv')
write.csv(sig_results, 'results/significant_circRNAs.csv')
