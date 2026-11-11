# GSE192849: Reproducible volcano plot

res <- read.csv(
  'results/all_differential_expression_results.csv',
  row.names = 1
)

# Retain finite, positive adjusted p-values
valid <- is.finite(res$log2FoldChange) &
  is.finite(res$padj) & res$padj > 0
d <- res[valid, , drop = FALSE]

category <- ifelse(
  d$padj < 0.05 & d$log2FoldChange > 0,
  'Upregulated',
  ifelse(d$padj < 0.05 & d$log2FoldChange < 0,
         'Downregulated', 'Not significant')
)

point_colors <- c(
  'Upregulated' = 'red',
  'Downregulated' = 'blue',
  'Not significant' = 'grey70'
)

dir.create('plots', showWarnings = FALSE)

png('plots/Volcano_plot_reproducible.png',
    width = 1200, height = 900, res = 150)

plot(
  d$log2FoldChange, -log10(d$padj),
  pch = 16, cex = 0.5,
  col = point_colors[category],
  xlab = 'Log2 fold change (Tumor vs Normal)',
  ylab = '-log10 adjusted p-value',
  main = 'Breast Cancer circRNA: Volcano Plot'
)

abline(h = -log10(0.05), lty = 2)
legend('topright',
       legend = names(point_colors),
       col = point_colors, pch = 16, cex = 0.8)
dev.off()

cat('Plottable circRNAs:', nrow(d), '\n')
cat('Significant circRNAs:',
    sum(d$padj < 0.05), '\n')
