#!/usr/bin/env Rscript
# Wilcoxon baseline (D_HARNESS_DESIGN.md v2 section 11).
# raw -> TSS relative abundance; rarefied / boots -> rarefied counts.
# Two-sided asymptotic wilcox.test with continuity correction, BH within the
# intention-to-test set, q < 0.05.

source(file.path(Sys.getenv("HARNESS_R_DIR", "."), "lib_common.R"))
WORKDIR <- env_require("WORKDIR")
CELL_ID <- env_require("CELL_ID"); SCENARIO <- env_require("SCENARIO")
REP <- as.integer(env_num("REP")); B <- as.integer(env_num("B"))
BASE_SEED <- env_num("BASE_SEED"); CELL_INDEX <- as.integer(env_num("CELL_INDEX"))
ds <- dataset_seed(BASE_SEED, CELL_INDEX, REP)

meta <- read.delim(file.path(WORKDIR, "meta.tsv"))
grp <- factor(meta$group, levels = c("G1", "G2"))
timings <- file.path(WORKDIR, "timings.tsv")

fit_one <- function(mat, tss) {
  if (tss) mat <- sweep(mat, 2, pmax(colSums(mat), 1L), "/")
  taxa <- rownames(mat)
  eff <- rowMeans(mat[, grp == "G2", drop = FALSE]) - rowMeans(mat[, grp == "G1", drop = FALSE])
  p <- vapply(seq_len(nrow(mat)), function(i) {
    if (all(mat[i, ] == mat[i, 1])) return(NA_real_)   # constant across samples
    suppressWarnings(wilcox.test(mat[i, ] ~ grp, exact = FALSE, correct = TRUE)$p.value)
  }, numeric(1))
  q <- rep(NA_real_, length(p)); q[!is.na(p)] <- p.adjust(p[!is.na(p)], "BH")
  data.frame(taxon = taxa, effect = eff, p = p, q = q,
             detected = as.integer(!is.na(q) & q < 0.05), stringsAsFactors = FALSE)
}
as_rows <- function(fitdf, condition) {
  out <- empty_results(fitdf$taxon, CELL_ID, SCENARIO, REP, "wilcoxon", condition)
  out$effect <- fitdf$effect; out$p <- fitdf$p; out$q <- fitdf$q
  out$detected <- fitdf$detected
  out$reason <- ifelse(is.na(fitdf$p), "constant", "ok")
  out
}

set.seed(ds + SEED_OFF_METHOD)
rows <- list()
t0 <- proc.time()[3]
raw <- read_counts(file.path(WORKDIR, "raw.tsv"))
rows$raw <- as_rows(fit_one(raw, tss = TRUE), "raw")
append_timing(timings, CELL_ID, REP, "wilcoxon", "raw", NA, proc.time()[3] - t0)

t0 <- proc.time()[3]
rar <- read_counts(file.path(WORKDIR, "rarefied.tsv"))
rows$rarefied <- as_rows(fit_one(rar, tss = FALSE), "rarefied")
append_timing(timings, CELL_ID, REP, "wilcoxon", "rarefied", NA, proc.time()[3] - t0)

taxa <- rownames(raw)
reps_list <- vector("list", B)
for (b in seq_len(B)) {
  t0 <- proc.time()[3]
  bm <- read_counts(file.path(WORKDIR, "boots", sprintf("b%03d.tsv", b)))
  f <- fit_one(bm, tss = FALSE)
  reps_list[[b]] <- f
  append_timing(timings, CELL_ID, REP, "wilcoxon", "boots", b, proc.time()[3] - t0)
}
bs <- summarize_boots(reps_list, taxa)
bo <- empty_results(taxa, CELL_ID, SCENARIO, REP, "wilcoxon", "boots")
bo$boots_B <- B; bo$boots_ok <- bs$boots_ok
bo$boots_detect_count <- bs$boots_detect_count
bo$boots_q_median <- bs$boots_q_median; bo$boots_q_q25 <- bs$boots_q_q25
bo$boots_q_q75 <- bs$boots_q_q75; bo$boots_effect_median <- bs$boots_effect_median
bo$boots_sign_agree <- bs$boots_sign_agree
bo$effect <- bs$boots_effect_median
bo$detected <- as.integer(bs$boots_ok > 0 & bs$boots_detect_count / pmax(bs$boots_ok, 1L) >= 0.5)
rows$boots <- bo

write_results(do.call(rbind, rows), file.path(WORKDIR, "results_wilcoxon.tsv"))
cat(sprintf("wilcoxon done: %s rep %d\n", CELL_ID, REP))
