#!/usr/bin/env Rscript
# B-3 (CG-08 部分採用): quantify the closure-induced relative shift of
# non-injected taxa and contrast it with the detection floor.
#
# Theory: injection multiplies spiked latent abundances by f (up) or 1/f
# (down); the total-mass factor c makes every NON-spiked taxon's relative
# abundance shrink by exactly 1/c (truth.tsv stores rel_fc_analytic = 1/c
# for spiked = 0). The estimand fixes these as absolute-null; B-3 asks
# whether the closure shift (log c) is anywhere near the effect sizes the
# methods actually flag (the "detection floor").
#
# Scales: ancombc2 effect = bias-corrected log fold change (natural log);
# maaslin3 effect = abundance-model coefficient (LOG transform = log2).
# The closure shift is converted into each method's own scale. wilcoxon's
# effect is a mean difference (not log scale) and is out of scope here.
# usage: Rscript b3_closure_shift.R <dsim_run_dir> <out_dir>

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
run_dir <- args[1]; out_dir <- args[2]
stopifnot(dir.exists(run_dir))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- Sys.glob(file.path(run_dir, "results", "spike_*", "rep_[0-9]*.tsv.gz"))
files <- files[!grepl("_truth", files)]
stopifnot(length(files) > 0)
message(sprintf("reading %d spike rep files ...", length(files)))

theory_rows <- list(); det_rows <- list()
for (f in files) {
  r <- read.delim(gzfile(f), stringsAsFactors = FALSE)
  r <- r[r$method %in% c("ancombc2", "maaslin3") &
         r$condition %in% c("raw", "rarefied"), ]
  t <- read.delim(gzfile(sub("\\.tsv\\.gz$", "_truth.tsv.gz", f)),
                  stringsAsFactors = FALSE)
  rel_null <- t$rel_fc_analytic[t$spiked == 0L][1]   # identical within rep
  theory_rows[[f]] <- data.frame(cell_id = r$cell_id[1], rep = r$rep[1],
                                 c_factor = 1 / rel_null)
  m <- merge(r, t[, c("taxon", "spiked", "tested")], by = "taxon", sort = FALSE)
  m <- m[m$tested == 1L & !is.na(m$detected) & m$detected == 1L &
         !is.na(m$effect), ]
  det_rows[[f]] <- m[, c("cell_id", "rep", "method", "condition",
                         "taxon", "effect", "spiked")]
}
theory <- do.call(rbind, theory_rows)
det <- do.call(rbind, det_rows)

th <- do.call(rbind, lapply(split(theory, theory$cell_id), function(g) data.frame(
  cell_id = g$cell_id[1], n_reps = nrow(g),
  c_mean = mean(g$c_factor), c_min = min(g$c_factor), c_max = max(g$c_factor),
  shift_ln_mean = mean(log(g$c_factor)),
  shift_log2_mean = mean(log2(g$c_factor)),
  stringsAsFactors = FALSE)))
rownames(th) <- NULL
write.table(th[order(th$cell_id), ],
            file.path(out_dir, "b3_closure_theory_by_cell.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

p05 <- function(x) quantile(x, 0.05, names = FALSE)
key <- c("cell_id", "method", "condition")
fl <- do.call(rbind, lapply(split(det, interaction(det[key], drop = TRUE)), function(g) {
  shift <- if (g$method[1] == "ancombc2") th$shift_ln_mean[th$cell_id == g$cell_id[1]]
           else th$shift_log2_mean[th$cell_id == g$cell_id[1]]
  fp <- abs(g$effect[g$spiked == 0L]); tp <- abs(g$effect[g$spiked == 1L])
  data.frame(g[1, key],
    effect_scale = ifelse(g$method[1] == "ancombc2", "ln", "log2"),
    closure_shift = shift,
    n_fp = length(fp), fp_abs_min = suppressWarnings(min(fp)),
    fp_abs_p05 = if (length(fp)) p05(fp) else NA_real_,
    fp_abs_median = if (length(fp)) median(fp) else NA_real_,
    frac_fp_below_shift = if (length(fp)) mean(fp <= shift) else NA_real_,
    floor_to_shift_ratio = if (length(fp)) p05(fp) / shift else NA_real_,
    n_tp = length(tp),
    tp_abs_median = if (length(tp)) median(tp) else NA_real_,
    stringsAsFactors = FALSE)
}))
rownames(fl) <- NULL
fl <- fl[order(fl$cell_id, fl$method, fl$condition), ]
write.table(fl, file.path(out_dir, "b3_detection_floor_vs_closure.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("theory cells: ", nrow(th), "; floor rows: ", nrow(fl))
print(th[, c("cell_id", "c_mean", "shift_log2_mean")])
print(fl[fl$cell_id == "spike_r100_f3",
         c("method", "condition", "closure_shift", "fp_abs_p05",
           "fp_abs_median", "frac_fp_below_shift", "floor_to_shift_ratio")])
