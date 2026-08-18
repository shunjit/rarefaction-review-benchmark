#!/usr/bin/env Rscript
# Aggregate a D-real run (D_HARNESS_DESIGN.md v2 section 5).
#
# D-real has no ground truth: metrics are DESCRIPTIVE distributions over the
# R resampled splits (median / IQR / min / max), never SEs or tests — the 80%
# subsamples overlap heavily (design section 5, review condition).
#
# Per (cell x method x condition x scenario):
#   n_detected per rep, detection fraction (n_detected / n_tested),
#   share of reps with >= 1 detection ("any-detection rate"),
#   plus the boots threshold variants (70/80/90%) as in the D-sim aggregator.
# Scenarios are paired within rep (same subsample, same tables):
#   real_depth = stratified median-depth split, real_perm = restricted
#   permutation -> paired per-rep differences (depth - perm) are summarised.
# Achieved depth ratios (geninfo) are joined into the per-rep table so the
# "detections vs achieved ratio" relation can be plotted downstream.
#
# usage: Rscript aggregate_real_run.R <run_dir>

args <- commandArgs(trailingOnly = TRUE)
run_dir <- args[1]
stopifnot(dir.exists(run_dir))
out_dir <- file.path(run_dir, "summary")
dir.create(out_dir, showWarnings = FALSE)

res_files <- Sys.glob(file.path(run_dir, "results", "*", "rep_[0-9]*.tsv.gz"))
res_files <- res_files[!grepl("_truth", res_files)]
stopifnot(length(res_files) > 0)

read_one <- function(f) {
  r <- read.delim(gzfile(f), stringsAsFactors = FALSE)
  tf <- sub("\\.tsv\\.gz$", "_truth.tsv.gz", f)
  t <- read.delim(gzfile(tf), stringsAsFactors = FALSE)
  merge(r, t[, c("taxon", "stratum", "tested")], by = "taxon", sort = FALSE)
}
message(sprintf("reading %d rep files ...", length(res_files)))
d <- do.call(rbind, lapply(res_files, read_one))
d <- d[d$tested == 1L, ]                       # intention-to-test universe
d$detected[is.na(d$detected)] <- 0L

# --- boots threshold variants (same construction as D-sim aggregator) --------
thr_rows <- function(dd, thr) {
  b <- dd[dd$condition == "boots" & !is.na(dd$boots_ok) & dd$boots_ok > 0, ]
  if (nrow(b) == 0L) return(NULL)                # B=0 run: no boots to threshold
  b$detected <- as.integer(b$boots_detect_count / b$boots_ok >= thr)
  b$condition <- sprintf("boots@%02d", round(100 * thr))
  b
}
base <- d[d$condition %in% c("raw", "rarefied", "boots"), ]
# B=0 runs (phase 1: raw/rarefied only) still emit boots rows with boots_ok=0;
# those are "condition not run", not zero detections — drop them.
base <- base[!(base$condition == "boots" & (is.na(base$boots_ok) | base$boots_ok == 0L)), ]
d_all <- rbind(base, do.call(rbind, lapply(c(0.7, 0.8, 0.9), function(t) thr_rows(d, t))))

# --- per-rep detection counts -------------------------------------------------
key <- c("cell_id", "method", "condition", "scenario")
grp <- interaction(d_all[c(key, "rep")], drop = TRUE)
per_rep <- do.call(rbind, lapply(split(d_all, grp), function(g) data.frame(
  g[1, key], rep = g$rep[1],
  n_tested = nrow(g), n_detected = sum(g$detected),
  frac_detected = sum(g$detected) / nrow(g),
  any_detected = as.integer(any(g$detected == 1L)),
  stringsAsFactors = FALSE)))
rownames(per_rep) <- NULL

# join achieved depth ratios from geninfo
gi_files <- Sys.glob(file.path(run_dir, "results", "*", "rep_[0-9]*_geninfo.tsv"))
gi <- do.call(rbind, lapply(gi_files, read.delim))
gi_key <- paste(gi$cell_id, gi$rep)
m <- match(paste(per_rep$cell_id, per_rep$rep), gi_key)
per_rep$achieved_ratio_depth <- gi$achieved_ratio_depth[m]
per_rep$achieved_ratio_perm <- gi$achieved_ratio_perm[m]
per_rep$raref_depth <- gi$raref_depth[m]
per_rep$n_sub <- gi$n_sub[m]
write.table(per_rep[order(per_rep$cell_id, per_rep$method, per_rep$condition,
                          per_rep$scenario, per_rep$rep), ],
            file.path(out_dir, "per_rep_detections.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- descriptive summaries per cell x method x condition x scenario ----------
q1 <- function(x) quantile(x, 0.25, names = FALSE)
q3 <- function(x) quantile(x, 0.75, names = FALSE)
sgrp <- interaction(per_rep[key], drop = TRUE)
summ <- do.call(rbind, lapply(split(per_rep, sgrp), function(g) data.frame(
  g[1, key], n_reps = nrow(g),
  tested_n_median = median(g$n_tested),
  det_median = median(g$n_detected), det_q25 = q1(g$n_detected),
  det_q75 = q3(g$n_detected), det_min = min(g$n_detected), det_max = max(g$n_detected),
  frac_median = median(g$frac_detected), frac_q25 = q1(g$frac_detected),
  frac_q75 = q3(g$frac_detected),
  any_detection_rate = mean(g$any_detected),
  stringsAsFactors = FALSE)))
rownames(summ) <- NULL
write.table(summ[order(summ$cell_id, summ$method, summ$condition, summ$scenario), ],
            file.path(out_dir, "metrics_real_by_cell.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- paired within-rep differences: depth split minus restricted permutation --
wide_key <- c("cell_id", "method", "condition", "rep")
dep <- per_rep[per_rep$scenario == "real_depth", ]
per <- per_rep[per_rep$scenario == "real_perm", ]
mm <- match(interaction(dep[wide_key], drop = FALSE), interaction(per[wide_key], drop = FALSE))
pair <- dep[!is.na(mm), wide_key]
pair$diff_detected <- dep$n_detected[!is.na(mm)] - per$n_detected[mm[!is.na(mm)]]
pgrp <- interaction(pair[c("cell_id", "method", "condition")], drop = TRUE)
pair_summ <- do.call(rbind, lapply(split(pair, pgrp), function(g) data.frame(
  g[1, c("cell_id", "method", "condition")], n_pairs = nrow(g),
  diff_median = median(g$diff_detected), diff_q25 = q1(g$diff_detected),
  diff_q75 = q3(g$diff_detected),
  frac_pairs_depth_higher = mean(g$diff_detected > 0),
  frac_pairs_equal = mean(g$diff_detected == 0),
  stringsAsFactors = FALSE)))
rownames(pair_summ) <- NULL
write.table(pair_summ[order(pair_summ$cell_id, pair_summ$method, pair_summ$condition), ],
            file.path(out_dir, "paired_depth_vs_perm.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- timings ------------------------------------------------------------------
t_files <- Sys.glob(file.path(run_dir, "results", "*", "rep_[0-9]*_timings.tsv"))
tt <- do.call(rbind, lapply(t_files, read.delim))
tg <- interaction(tt$method, tt$condition, drop = TRUE)
tsum <- do.call(rbind, lapply(split(tt, tg), function(g) data.frame(
  method = g$method[1], condition = g$condition[1], n_fits = nrow(g),
  sec_median = median(g$elapsed_sec), sec_mean = mean(g$elapsed_sec),
  sec_total = sum(g$elapsed_sec), stringsAsFactors = FALSE)))
tsum <- tsum[order(tsum$method, tsum$condition), ]
rownames(tsum) <- NULL
write.table(tsum, file.path(out_dir, "timings_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message(sprintf("aggregated %d reps x scenario rows over %d cells",
                length(unique(paste(d$cell_id, d$rep))), length(unique(d$cell_id))))
print(summ[summ$condition %in% c("raw", "rarefied", "boots"), ])
