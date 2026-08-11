# Clustering parameters are compared, establishing 50 percent identity at coverage mode 1 as primary.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/run_happi_focal.R
# Output: work/focal_genus_pangenome/matrices/happi_results_<clustering>.tsv
# RUN_HAPPI_FOCAL_V2_20260804
suppressPackageStartupMessages({ library(happi) })

ROOT <- "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
WORK <- file.path(ROOT, "work/focal_genus_pangenome/matrices")
TAG  <- Sys.getenv("CLUSTERING", "id80_cov1")
MAT  <- file.path(WORK, sprintf("presence_%s.tsv", TAG))
MET  <- file.path(WORK, "unit_metadata.tsv")
OUT  <- file.path(WORK, sprintf("happi_results_%s.tsv", TAG))

cat("clustering:", TAG, "\n")
stopifnot(file.exists(MAT), file.exists(MET))

pres <- read.delim(MAT, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.delim(MET, stringsAsFactors = FALSE)
genomes <- colnames(pres)[-1]
meta <- meta[match(genomes, meta$genome), ]
stopifnot(all(meta$genome == genomes))

meta$amphibian <- as.integer(meta$group == "amphibian")
meta$comp <- meta$completeness / 100

cat("clusters:", nrow(pres), " units:", nrow(meta), "\n")
print(table(meta$group))
cat("completeness: amphibian mean", round(mean(meta$comp[meta$amphibian == 1]) * 100, 1),
    " reference mean", round(mean(meta$comp[meta$amphibian == 0]) * 100, 1), "\n\n")

# V1 passed column-name strings; happi calls nrow(covariate) so it needs a
# design matrix and a numeric vector. pbLRT is unavailable here (it requires
# arguments the fit object does not retain), so the asymptotic pvalue_LRT is
# used. That is less well calibrated than the parametric bootstrap the happi
# paper prefers for small samples, and must be stated as a limitation.
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
    cluster = pres$cluster[i],
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
cat("\nfitted", nfit, " skipped (invariant)", nskip, " failed", nfail, "\n")
res$q <- p.adjust(res$p, method = "BH")
res <- res[order(res$p), ]
write.table(res, OUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", OUT, "\n\n")

sig <- !is.na(res$q) & res$q < 0.05
cat("BH q < 0.05:", sum(sig), "\n")
cat("  higher prevalence in amphibian:", sum(sig & res$diff > 0), "\n")
cat("  higher prevalence in reference:", sum(sig & res$diff < 0), "\n\n")
cat("top 15 by p:\n")
print(head(res[, c("cluster", "prev_amphibian", "prev_reference", "diff", "p", "q")], 15),
      row.names = FALSE)
cat("\nhappi is conservative when genome quality correlates with the group\n")
cat("being tested, which it does here. Null results are weak evidence of no\n")
cat("difference rather than evidence of no difference.\n")
# SENTINEL_END
