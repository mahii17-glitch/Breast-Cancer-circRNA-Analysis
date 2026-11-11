# GSE192849: Load and validate circRNA counts

counts_raw <- read.delim(
  'data/GSE192849_raw_counts.txt.gz',
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sample_names <- c(
  'Tumor1', 'Normal1',
  'Tumor2', 'Normal2',
  'Tumor3', 'Normal3'
)

stopifnot(
  'circRNA' %in% names(counts_raw),
  all(sample_names %in% names(counts_raw))
)

counts_int <- as.matrix(
  counts_raw[, sample_names, drop = FALSE]
)

rownames(counts_int) <- counts_raw$circRNA

stopifnot(
  nrow(counts_int) == 47433,
  !anyNA(counts_int),
  !anyDuplicated(rownames(counts_int)),
  is.numeric(counts_int),
  all(counts_int >= 0),
  all(counts_int == floor(counts_int))
)

storage.mode(counts_int) <- 'integer'

cat('circRNAs:', nrow(counts_int), '\n')
cat('Samples:', ncol(counts_int), '\n')
cat('Validation: PASSED\n')
