# Reproducible circRNA genomic containment analysis
# Uses saved gene coordinates; no internet required

circ <- read.csv(
  'results/circRNA_genomic_coordinates.csv',
  check.names = FALSE
)

genes <- read.csv(
  'data/annotation/GRCh38_gene_coordinates.csv',
  check.names = FALSE
)

matches <- lapply(seq_len(nrow(circ)), function(i) {
  x <- circ[i, ]
  hits <- genes[
    genes$chromosome == x$chromosome &
    genes$start <= x$start &
    genes$end >= x$end &
    genes$strand == x$strand,
    c('gene_id', 'gene_name'),
    drop = FALSE
  ]
  data.frame(circRNA_ID = x$circRNA_ID, hits)
})

verified <- do.call(rbind, matches)
rownames(verified) <- NULL

dir.create('results', showWarnings = FALSE)
write.csv(
  verified,
  'results/verified_candidate_host_genes.csv',
  row.names = FALSE
)

cat('Verified annotation records:',
    nrow(verified), '\n')

# Note: genome assembly and coordinate convention
# of the original circRNA dataset require confirmation.
# Genomic containment does not prove host-gene origin.
