# Orthologous-group prevalence modelled against completeness for the amphibian clade against the other Angelakisella genomes.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/run_happi_angelakisella.R
# Output: results/happi_angelakisella_og.tsv and permutation outputs
# RUN_HAPPI_ANGELAKISELLA_V1
# Clade-content test for the Angelakisella neighborhood: amphibian clade
# (31) vs non-amphibian Angelakisella (26), eggNOG Bacteria-level OG
# presence, completeness-aware via happi, permutation null via env vars.
# Cloned from run_happi_merged.R; happi call parameters unchanged.
suppressPackageStartupMessages({ library(happi) })
ROOT <- "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
MAT  <- file.path(ROOT, "results/angelakisella_matrix_og_v2.tsv")
MET  <- file.path(ROOT, "results/angelakisella_happi_metadata.tsv")
PERM <- as.integer(Sys.getenv("PERMUTE", "0"))
SEED <- as.integer(Sys.getenv("PERM_SEED", "0"))
OUT <- if (PERM == 1)
  file.path(ROOT, sprintf("results/happi_angelakisella_og_perm%d.tsv", SEED)) else
  file.path(ROOT, "results/happi_angelakisella_og.tsv")
if (file.exists(OUT)) stop("STOP: output exists, refusing to overwrite: ", OUT)
stopifnot(file.exists(MAT), file.exists(MET))
pres <- read.delim(MAT, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.delim(MET, stringsAsFactors = FALSE)
keep <- c("unit", meta$genome)
stopifnot(all(meta$genome %in% colnames(pres)))
pres <- pres[, keep]
stopifnot(all(colnames(pres)[-1] == meta$genome))
meta$amphibian <- as.integer(meta$group == "amphibian")
meta$comp <- meta$completeness / 100
if (PERM == 1) {
  set.seed(SEED)
  meta$amphibian <- sample(meta$amphibian)
  cat("LABELS PERMUTED. This run is a null control.\n")
}
cat("units:", nrow(pres), " genomes:", nrow(meta), "\n")
print(table(meta$group))
cat("completeness: amphibian",
    round(mean(meta$comp[meta$amphibian == 1]) * 100, 1),
    " non_amphibian",
    round(mean(meta$comp[meta$amphibian == 0]) * 100, 1), "\n\n")
X <- model.matrix(~ amphibian, data = meta)
amph <- meta$amphibian == 1
res <- vector("list", nrow(pres))
nfit <- 0; nskip <- 0; nfail <- 0
t0 <- Sys.time()
for (i in seq_len(nrow(pres))) {
  y <- as.integer(pres[i, -1])
  if (length(unique(y)) < 2) { nskip <- nskip + 1; next }
  f <- try(happi(outcome = y, covariate = X, quality_var = meta$comp,
                 method = "splines", firth = TRUE, spline_df = 3,
                 max_iterations = 500, change_threshold = 0.05,
                 epsilon = 0, verbose = FALSE), silent = TRUE)
  if (inherits(f, "try-error") || is.null(f$pvalue_LRT)) { nfail <- nfail + 1; next }
  nfit <- nfit + 1
  res[[i]] <- data.frame(
    og = pres$unit[i],
    n_amphibian = sum(y[amph]), n_non_amphibian = sum(y[!amph]),
    prev_amphibian = mean(y[amph]), prev_non_amphibian = mean(y[!amph]),
    diff = mean(y[amph]) - mean(y[!amph]),
    pvalue = f$pvalue_LRT)
  if (i %% 250 == 0)
    cat(i, "of", nrow(pres), " fitted:", nfit, " elapsed:",
        round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
}
res <- do.call(rbind, res[!sapply(res, is.null)])
res$qvalue <- p.adjust(res$pvalue, method = "BH")
res <- res[order(res$pvalue), ]
write.table(res, OUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nfitted:", nfit, " skipped invariant:", nskip, " failed:", nfail, "\n")
cat("significant at BH q<0.05:", sum(res$qvalue < 0.05), "\n")
cat("wrote:", OUT, "\n")
cat("HAPPI_ANGELAKISELLA_FINISHED perm=", PERM, " seed=", SEED, "\n")
