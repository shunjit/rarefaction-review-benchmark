#!/usr/bin/env Rscript
# B=100 convergence check (D_HARNESS_DESIGN.md v2 section 3).
#
# The main run reduced boots to B=20 (main cells) / B=10 (sensitivity cells)
# via the pre-registered fallback; section 3 requires confirming at B=100, on
# one null and one spike-in cell, that detection fractions, ranks and binary
# calls are converged. The convergence run (cells_boots_convergence.tsv,
# cell_index 40/41) draws INDEPENDENT datasets (different cell_index => a
# different seed stream), so all comparisons are at the distribution / metric
# level — never per-taxon pairing across runs.
#
# Three registers compared between B=20 (main run cells, 150 reps) and B=100
# (conv run cells, 50 reps), same code path for both:
#   1. detection fraction: distribution of the per-taxon stability score
#      (boots_detect_count / boots_ok), plus the mass of taxa whose score sits
#      within the binomial 95% band around the 0.5 call threshold — the taxa
#      whose binary call can flip due to B alone
#   2. rank: per-rep AUC of the score separating spiked from null taxa
#      (spike pair only; Mann-Whitney statistic / (n1*n0))
#   3. binary calls: cell-level metrics of the thresholded call at
#      50/70/80/90% — null: taxon-FPR; spike: mean FDP + power_cot
#
# usage: Rscript convergence_b100.R <main_run_dir> <conv_run_dir> <out_dir>

args <- commandArgs(trailingOnly = TRUE)
main_dir <- args[1]; conv_dir <- args[2]; out_dir <- args[3]
stopifnot(dir.exists(main_dir), dir.exists(conv_dir))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

PAIRS <- data.frame(
  pair = c("null_r100", "spike_r100_f3"),
  main_cell = c("null_r100", "spike_r100_f3"),
  conv_cell = c("conv_null_r100_B100", "conv_spike_r100_f3_B100"),
  stringsAsFactors = FALSE)
METHODS <- c("ancombc2", "maaslin3", "wilcoxon")

read_boots <- function(run_dir, cell) {
  fs <- Sys.glob(file.path(run_dir, "results", cell, "rep_[0-9]*.tsv.gz"))
  fs <- fs[!grepl("_truth", fs)]
  stopifnot(length(fs) > 0)
  do.call(rbind, lapply(fs, function(f) {
    r <- read.delim(gzfile(f), stringsAsFactors = FALSE)
    t <- read.delim(gzfile(sub("\\.tsv\\.gz$", "_truth.tsv.gz", f)),
                    stringsAsFactors = FALSE)
    r <- r[r$condition == "boots" & r$method %in% METHODS &
             !is.na(r$boots_ok) & r$boots_ok > 0, ]
    r <- merge(r, t[, c("taxon", "spiked", "tested")], by = "taxon", sort = FALSE)
    r <- r[r$tested == 1L, ]
    data.frame(rep = r$rep, method = r$method, taxon = r$taxon,
               spiked = r$spiked, score = r$boots_detect_count / r$boots_ok,
               B_ok = r$boots_ok, stringsAsFactors = FALSE)
  }))
}

q1 <- function(x) quantile(x, 0.25, names = FALSE)
q3 <- function(x) quantile(x, 0.75, names = FALSE)

summ_rows <- list(); auc_rows <- list(); call_rows <- list()
for (i in seq_len(nrow(PAIRS))) {
  p <- PAIRS[i, ]
  for (src in c("B20", "B100")) {
    d <- if (src == "B20") read_boots(main_dir, p$main_cell)
         else read_boots(conv_dir, p$conv_cell)
    is_spike <- grepl("spike", p$pair)
    for (m in METHODS) {
      g <- d[d$method == m, ]
      B_nom <- median(g$B_ok)
      band <- 1.96 * sqrt(0.25 / B_nom)          # binomial SE at p=0.5
      n_reps <- length(unique(g$rep))
      # -- register 1: score distribution + call-flippable mass ---------------
      summ_rows[[length(summ_rows) + 1]] <- data.frame(
        pair = p$pair, method = m, source = src, n_reps = n_reps,
        B_nominal = B_nom,
        score_pos_frac = mean(g$score > 0),
        score_median_pos = if (any(g$score > 0)) median(g$score[g$score > 0]) else NA_real_,
        score_q3_pos = if (any(g$score > 0)) q3(g$score[g$score > 0]) else NA_real_,
        flip_band_halfwidth = band,
        flippable_mass = mean(abs(g$score - 0.5) < band),
        stringsAsFactors = FALSE)
      # -- register 2: rank (spike only): per-rep AUC spiked vs null ----------
      if (is_spike) {
        aucs <- vapply(split(g, g$rep), function(r) {
          x <- r$score[r$spiked == 1L]; y <- r$score[r$spiked == 0L]
          if (!length(x) || !length(y)) return(NA_real_)
          rk <- rank(c(x, y))
          (sum(rk[seq_along(x)]) - length(x) * (length(x) + 1) / 2) /
            (length(x) * length(y))
        }, numeric(1))
        auc_rows[[length(auc_rows) + 1]] <- data.frame(
          pair = p$pair, method = m, source = src, n_reps = n_reps,
          auc_median = median(aucs, na.rm = TRUE),
          auc_q25 = q1(aucs[!is.na(aucs)]), auc_q75 = q3(aucs[!is.na(aucs)]),
          stringsAsFactors = FALSE)
      }
      # -- register 3: binary calls at 50/70/80/90% ---------------------------
      for (thr in c(0.5, 0.7, 0.8, 0.9)) {
        g$call <- as.integer(g$score >= thr)
        per_rep <- split(g, g$rep)
        if (is_spike) {
          fdp <- vapply(per_rep, function(r)
            sum(r$call & r$spiked == 0L) / max(sum(r$call), 1), numeric(1))
          pow <- vapply(per_rep, function(r)
            sum(r$call & r$spiked == 1L) / max(sum(r$spiked == 1L), 1), numeric(1))
          call_rows[[length(call_rows) + 1]] <- data.frame(
            pair = p$pair, method = m, source = src, thr = thr,
            fpr_taxon = sum(g$call & g$spiked == 0L) / sum(g$spiked == 0L),
            mean_fdp = mean(fdp), power_cot = mean(pow), stringsAsFactors = FALSE)
        } else {
          call_rows[[length(call_rows) + 1]] <- data.frame(
            pair = p$pair, method = m, source = src, thr = thr,
            fpr_taxon = sum(g$call) / nrow(g),
            mean_fdp = NA_real_, power_cot = NA_real_, stringsAsFactors = FALSE)
        }
      }
    }
  }
}

s <- do.call(rbind, summ_rows); a <- do.call(rbind, auc_rows)
cc <- do.call(rbind, call_rows)
write.table(s, file.path(out_dir, "convergence_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(a, file.path(out_dir, "convergence_auc.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cc[order(cc$pair, cc$method, cc$thr, cc$source), ],
            file.path(out_dir, "convergence_calls.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("=== binary-call metrics, B20 vs B100 (key register) ===")
print(cc[order(cc$pair, cc$method, cc$thr, cc$source), ], row.names = FALSE)
message("=== rank AUC (spike) ===")
print(a, row.names = FALSE)
message("=== flippable mass around 0.5 ===")
print(s[, c("pair", "method", "source", "B_nominal", "flippable_mass")],
      row.names = FALSE)
