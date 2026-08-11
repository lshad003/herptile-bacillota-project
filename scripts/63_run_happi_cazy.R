# Detectable prevalence is modelled against completeness for the CAZy families.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/run_happi_cazy_v2.R
# Output: results/happi_cazy_focal.tsv and permuted controls
# RUN_HAPPI_CAZY_V2_20260806
suppressPackageStartupMessages({ library(happi) })

ROOT <- "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
MTX  <- file.path(ROOT, "results/cazy_focal_family_matrix.tsv")
QUAL <- file.path(ROOT, "work/focal_genus_pangenome/checkm2_out/quality_report.tsv")
PERM <- as.integer(Sys.getenv("PERMUTE", "0"))
SEED <- as.integer(Sys.getenv("PERM_SEED", "0"))
OUT  <- if (PERM == 1)
  file.path(ROOT, sprintf("results/happi_cazy_focal_perm%d.tsv", SEED)) else
  file.path(ROOT, "results/happi_cazy_focal.tsv")
MIN_G <- 10

stopifnot(file.exists(MTX), file.exists(QUAL))
if (file.exists(OUT)) stop("refusing to overwrite ", OUT)
cat("permuted:", PERM, " seed:", SEED, "\n\n")

m <- read.delim(MTX, check.names = FALSE, stringsAsFactors = FALSE)
m <- m[m$set == "focal125" & m$arm %in% c("amphibian", "reference"), ]
cat("genomes:", nrow(m), " (endotherm-labelled genome excluded)\n")
print(table(m$genus, m$arm))

q <- read.delim(QUAL, check.names = FALSE, stringsAsFactors = FALSE)
gid <- vapply(strsplit(sub("\\.faa$", "", basename(q$Name)), "__", fixed = TRUE),
              function(z) z[length(z)], character(1))
qi <- match(m$genome_id, gid)
stopifnot(!any(is.na(qi)))
m$comp <- as.numeric(q$Completeness[qi]) / 100
cat("completeness: amphibian", round(mean(m$comp[m$arm == "amphibian"]) * 100, 1),
    " reference", round(mean(m$comp[m$arm == "reference"]) * 100, 1), "\n\n")

m$amphibian <- as.integer(m$arm == "amphibian")
if (PERM == 1) {
  set.seed(SEED)
  m$amphibian <- sample(m$amphibian)
  cat("LABELS PERMUTED. This run is a null control.\n\n")
}

meta_cols <- c("genome_id", "set", "genus", "arm", "clade", "total", "comp", "amphibian")
fams <- setdiff(names(m), meta_cols)
pres <- sapply(fams, function(f) sum(m[[f]] > 0))
keep <- fams[pres >= MIN_G & pres <= (nrow(m) - MIN_G)]
cat("families:", length(fams), " testable:", length(keep), "\n")

X <- model.matrix(~ amphibian + genus, data = m)
cat("design:", paste(colnames(X), collapse = ", "), "\n\n")
amph <- m$amphibian == 1

res <- vector("list", length(keep))
nfit <- 0; nfail <- 0; firsterr <- ""
t0 <- Sys.time()
for (i in seq_along(keep)) {
  y <- as.integer(m[[keep[i]]] > 0)
  if (length(unique(y)) < 2) { nfail <- nfail + 1; next }
  f <- try(happi(outcome = y, covariate = X, quality_var = m$comp,
                 method = "splines", firth = TRUE, spline_df = 3,
                 max_iterations = 500, change_threshold = 0.05,
                 epsilon = 0, verbose = FALSE), silent = TRUE)
  if (inherits(f, "try-error")) {
    nfail <- nfail + 1
    if (firsterr == "") firsterr <- as.character(f)
    next
  }
  if (is.null(f$pvalue_LRT)) { nfail <- nfail + 1; next }
  nfit <- nfit + 1
  res[[i]] <- data.frame(
    family = keep[i],
    n_amphibian = sum(y[amph]), n_reference = sum(y[!amph]),
    prev_amphibian = mean(y[amph]), prev_reference = mean(y[!amph]),
    diff = mean(y[amph]) - mean(y[!amph]),
    beta = if (length(f$beta) >= 2) f$beta[2] else NA_real_,
    LRT = f$LRT, p = f$pvalue_LRT, stringsAsFactors = FALSE)
  if (nfit %% 50 == 0)
    cat(sprintf("  %d fitted, %.1f min\n", nfit,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
if (nfit == 0) {
  cat("\nFIRST ERROR MESSAGE:\n"); cat(firsterr, "\n")
  stop("nothing fitted")
}
res <- do.call(rbind, res)
cat("\nfitted", nfit, " failed", nfail, "\n")
if (nfail > 0 && firsterr != "") { cat("first error:\n"); cat(firsterr, "\n") }

res$q <- p.adjust(res$p, method = "BH")
res <- res[order(res$p), ]
write.table(res, OUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", OUT, "\n\n")

sig <- !is.na(res$q) & res$q < 0.05
if (PERM == 1) {
  cat(sprintf("PERMUTED: BH q<0.05: %d of %d (%.1f%%)\n",
              sum(sig), nrow(res), 100 * sum(sig) / nrow(res)))
  cat("ANY NON-ZERO COUNT HERE INVALIDATES THE REAL RUN.\n")
} else {
  cat("BH q < 0.05:", sum(sig), "of", nrow(res), "\n")
  cat("  higher in amphibian:", sum(sig & res$diff > 0), "\n")
  cat("  higher in reference:", sum(sig & res$diff < 0), "\n\n")
  cat("top 25 by p:\n")
  print(head(res[, c("family", "prev_amphibian", "prev_reference", "diff", "p", "q")], 25),
        row.names = FALSE, digits = 3)
  cat("\nGenus is adjusted for, not tested. happi is conservative when\n")
  cat("quality correlates with group, which it does here.\n")
}
# RUN_HAPPI_CAZY_V2_20260806_END
