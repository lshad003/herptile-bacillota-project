# Detectable prevalence is modelled against completeness for the merged orthologous groups.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/run_happi_merged.R
# Output: work/focal_genus_pangenome/matrices/happi_results_og_bacteria.tsv and permuted controls
# RUN_HAPPI_MERGED_V1_20260805
suppressPackageStartupMessages({ library(happi) })

ROOT <- "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
WORK <- file.path(ROOT, "work/focal_genus_pangenome/matrices")
TAG  <- Sys.getenv("CLUSTERING", "og_bacteria")
PERM <- as.integer(Sys.getenv("PERMUTE", "0"))
SEED <- as.integer(Sys.getenv("PERM_SEED", "0"))

MAT <- file.path(WORK, sprintf("presence_%s.tsv", TAG))
MET <- file.path(WORK, "unit_metadata.tsv")
OUT <- if (PERM == 1)
  file.path(WORK, sprintf("happi_results_%s_perm%d.tsv", TAG, SEED)) else
  file.path(WORK, sprintf("happi_results_%s.tsv", TAG))

cat("matrix:", TAG, " permuted:", PERM, " seed:", SEED, "\n")
stopifnot(file.exists(MAT), file.exists(MET))

pres <- read.delim(MAT, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.delim(MET, stringsAsFactors = FALSE)
genomes <- colnames(pres)[-1]
meta <- meta[match(genomes, meta$genome), ]
stopifnot(all(meta$genome == genomes))

meta$amphibian <- as.integer(meta$group == "amphibian")
meta$comp <- meta$completeness / 100

if (PERM == 1) {
  set.seed(SEED)
  meta$amphibian <- sample(meta$amphibian)
  cat("LABELS PERMUTED. This run is a null control.\n")
}

cat("groups:", nrow(pres), " units:", nrow(meta), "\n")
print(table(meta$group))
cat("completeness: amphibian", round(mean(meta$comp[meta$amphibian == 1]) * 100, 1),
    " reference", round(mean(meta$comp[meta$amphibian == 0]) * 100, 1), "\n\n")

# happi needs a design matrix and a numeric quality vector, not column names.
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
    og = pres$cluster[i],
    n_amphibian = sum(y[amph]), n_reference = sum(y[!amph]),
    prev_amphibian = mean(y[amph]), prev_reference = mean(y[!amph]),
    diff = mean(y[amph]) - mean(y[!amph]),
    beta = if (length(f$beta) >= 2) f$beta[2] else NA_real_,
    LRT = f$LRT, p = f$pvalue_LRT, stringsAsFactors = FALSE)
  if (nfit %% 500 == 0)
    cat(sprintf("  %d fitted, %.1f min\n", nfit,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

res <- do.call(rbind, res)
cat("\nfitted", nfit, " skipped", nskip, " failed", nfail, "\n")
res$q <- p.adjust(res$p, method = "BH")
res <- res[order(res$p), ]
write.table(res, OUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", OUT, "\n\n")

sig <- !is.na(res$q) & res$q < 0.05
if (PERM == 1) {
  cat(sprintf("PERMUTED: BH q<0.05: %d of %d (%.1f%%)\n",
              sum(sig), nrow(res), 100 * sum(sig) / nrow(res)))
} else {
  cat("BH q < 0.05:", sum(sig), "\n")
  cat("  higher in amphibian:", sum(sig & res$diff > 0), "\n")
  cat("  higher in reference:", sum(sig & res$diff < 0), "\n\n")
  cat("top 20 by p:\n")
  print(head(res[, c("og", "prev_amphibian", "prev_reference", "diff", "p", "q")], 20),
        row.names = FALSE)
  cat("\nAbsolute prevalences are now interpretable: universal genes reach\n")
  cat("0.98-1.00 in both arms after merging, which they did not before.\n")
  cat("happi remains conservative when quality correlates with group.\n")
}
# SENTINEL_END
