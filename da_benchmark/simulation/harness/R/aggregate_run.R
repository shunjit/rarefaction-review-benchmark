#!/usr/bin/env Rscript
# Aggregate a D-sim run: metrics per cell x method x condition (+ boots
# threshold sensitivity) and a timing summary for 判断2.
# Metrics follow D_HARNESS_DESIGN.md v2 section 10:
#   null : taxon-level FPR, mean FP/dataset, FWER, mean FDP
#   spike: empirical FDR (+MC SE), power (intention-to-test and
#          conditional-on-tested), per-stratum FDR/power
# usage: Rscript aggregate_run.R <run_dir>

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
  merge(r, t[, c("taxon", "stratum", "spiked", "tested")], by = "taxon", sort = FALSE)
}
message(sprintf("reading %d rep files ...", length(res_files)))
d <- do.call(rbind, lapply(res_files, read_one))
d <- d[d$tested == 1L, ]                     # intention-to-test universe
d$detected[is.na(d$detected)] <- 0L
d$is_null <- d$spiked == 0L

# --- boots threshold variants ------------------------------------------------
thr_rows <- function(dd, thr) {
  b <- dd[dd$condition == "boots" & dd$boots_ok > 0, ]
  b$detected <- as.integer(b$boots_detect_count / b$boots_ok >= thr)
  b$condition <- sprintf("boots@%02d", round(100 * thr))
  b
}
base <- d[d$condition %in% c("raw", "rarefied", "boots"), ]
d_all <- rbind(base, do.call(rbind, lapply(c(0.7, 0.8, 0.9), function(t) thr_rows(d, t))))

key <- c("cell_id", "scenario", "method", "condition")
grp <- interaction(d_all[key], drop = TRUE)
summ <- do.call(rbind, lapply(split(d_all, grp), function(g) {
  per_rep <- split(g, g$rep)
  fdp <- vapply(per_rep, function(r) {
    disc <- sum(r$detected); fp <- sum(r$detected & r$is_null); fp / max(disc, 1)
  }, numeric(1))
  fp_n <- vapply(per_rep, function(r) sum(r$detected & r$is_null), numeric(1))
  any_fp <- vapply(per_rep, function(r) any(r$detected & r$is_null), logical(1))
  row <- data.frame(g[1, key],
    n_reps = length(per_rep),
    tested_n_mean = nrow(g) / length(per_rep),
    fpr_taxon = sum(g$detected & g$is_null) / sum(g$is_null),
    mean_fp = mean(fp_n), fwer = mean(any_fp),
    mean_fdp = mean(fdp), fdr_mc_se = sd(fdp) / sqrt(length(per_rep)),
    power_itt = NA_real_, power_cot = NA_real_,
    stringsAsFactors = FALSE)
  if (g$scenario[1] == "spike") {
    tp <- vapply(per_rep, function(r) sum(r$detected & !r$is_null), numeric(1))
    spiked_tested <- vapply(per_rep, function(r) sum(!r$is_null), numeric(1))
    spiked_all <- vapply(per_rep, function(r) round(0.10 * 300), numeric(1))  # placeholder, overwritten below
    row$power_cot <- mean(tp / pmax(spiked_tested, 1))
    row$power_itt <- NA_real_   # filled from geninfo below
  }
  row
}))
rownames(summ) <- NULL

# intention-to-test power denominators come from geninfo (spiked_n per rep)
gi_files <- Sys.glob(file.path(run_dir, "results", "*", "rep_[0-9]*_geninfo.tsv"))
gi <- do.call(rbind, lapply(gi_files, read.delim))
if (any(gi$scenario == "spike")) {
  sp <- d_all[d_all$scenario == "spike", ]
  gsp <- interaction(sp[key], drop = TRUE)
  itt <- vapply(split(sp, gsp), function(g) {
    per_rep <- split(g, g$rep)
    tp <- vapply(per_rep, function(r) sum(r$detected & !r$is_null), numeric(1))
    denom <- gi$spiked_n[match(paste(g$cell_id[1], names(per_rep)),
                               paste(gi$cell_id, gi$rep))]
    mean(tp / denom)
  }, numeric(1))
  summ$power_itt[match(names(itt), interaction(summ[key], drop = TRUE))] <- itt
}
write.table(summ[order(summ$cell_id, summ$method, summ$condition), ],
            file.path(out_dir, "metrics_by_cell.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- per-stratum metrics (spike cells) --------------------------------------
sp <- d_all[d_all$scenario == "spike", ]
if (nrow(sp)) {
  gs <- interaction(sp[c(key, "stratum")], drop = TRUE)
  strat <- do.call(rbind, lapply(split(sp, gs), function(g) {
    per_rep <- split(g, g$rep)
    fdp <- vapply(per_rep, function(r) sum(r$detected & r$is_null) / max(sum(r$detected), 1), numeric(1))
    pow <- vapply(per_rep, function(r) sum(r$detected & !r$is_null) / max(sum(!r$is_null), 1), numeric(1))
    data.frame(g[1, c(key, "stratum")], n_reps = length(per_rep),
               fdr = mean(fdp), power_cot = mean(pow), stringsAsFactors = FALSE)
  }))
  rownames(strat) <- NULL
  write.table(strat, file.path(out_dir, "metrics_by_stratum.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# --- timings -----------------------------------------------------------------
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

per_rep_cost <- sum(tt$elapsed_sec) / length(unique(paste(tt$cell_id, tt$rep)))
message(sprintf("aggregated %d reps; mean wall cost per rep (sum of fit+gen timings) = %.1f s",
                length(unique(paste(d$cell_id, d$rep))), per_rep_cost))
print(tsum)
