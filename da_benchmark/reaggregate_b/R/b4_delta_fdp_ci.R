#!/usr/bin/env Rscript
# B-4 (CG-09): paired per-rep bootstrap CI for the rarefied - raw difference
# in FDP (and, as the mandated power companion, in intention-to-test power),
# for the 9 spike cells x {ancombc2, maaslin3}. Pairing is exact: within a
# rep both conditions share the same generated dataset (設計書 §12).
# Default detection rule (q < 0.05). 10,000 bootstrap resamples over reps,
# percentile 95% CI, fixed seed 20260816.
# usage: Rscript b4_delta_fdp_ci.R <dsim_run_dir> <out_dir>

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
run_dir <- args[1]; out_dir <- args[2]
stopifnot(dir.exists(run_dir))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- Sys.glob(file.path(run_dir, "results", "spike_*", "rep_[0-9]*.tsv.gz"))
files <- files[!grepl("_truth", files)]
stopifnot(length(files) > 0)
message(sprintf("reading %d spike rep files ...", length(files)))

per_rep <- list()
for (f in files) {
  r <- read.delim(gzfile(f), stringsAsFactors = FALSE)
  r <- r[r$method %in% c("ancombc2", "maaslin3") &
         r$condition %in% c("raw", "rarefied"), ]
  t <- read.delim(gzfile(sub("\\.tsv\\.gz$", "_truth.tsv.gz", f)),
                  stringsAsFactors = FALSE)
  m <- merge(r, t[, c("taxon", "spiked", "tested")], by = "taxon", sort = FALSE)
  m <- m[m$tested == 1L, ]
  m$detected[is.na(m$detected)] <- 0L
  gi <- read.delim(sub("\\.tsv\\.gz$", "_geninfo.tsv", f))
  per_rep[[f]] <- do.call(rbind, lapply(
    split(m, interaction(m$method, m$condition, drop = TRUE)), function(g) {
      det <- g$detected == 1L
      data.frame(cell_id = g$cell_id[1], rep = g$rep[1], method = g$method[1],
        condition = g$condition[1],
        fdp = sum(det & g$spiked == 0L) / max(sum(det), 1),
        power_itt = sum(det & g$spiked == 1L) / gi$spiked_n[1],
        stringsAsFactors = FALSE)
    }))
}
pr <- do.call(rbind, per_rep)
rownames(pr) <- NULL
write.table(pr[order(pr$cell_id, pr$method, pr$condition, pr$rep), ],
            file.path(out_dir, "b4_per_rep_fdp_power.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

BOOT <- 10000L
set.seed(20260816)
boot_ci <- function(d) {
  n <- length(d)
  bm <- vapply(seq_len(BOOT), function(i) mean(d[sample.int(n, n, replace = TRUE)]),
               numeric(1))
  quantile(bm, c(0.025, 0.975), names = FALSE)
}

cells <- unique(pr$cell_id)
out <- list()
for (cid in cells) for (meth in c("ancombc2", "maaslin3")) {
  raw <- pr[pr$cell_id == cid & pr$method == meth & pr$condition == "raw", ]
  rar <- pr[pr$cell_id == cid & pr$method == meth & pr$condition == "rarefied", ]
  i <- match(raw$rep, rar$rep)
  stopifnot(!anyNA(i))
  d_fdp <- rar$fdp[i] - raw$fdp
  d_pow <- rar$power_itt[i] - raw$power_itt
  ci_f <- boot_ci(d_fdp); ci_p <- boot_ci(d_pow)
  out[[paste(cid, meth)]] <- data.frame(cell_id = cid, method = meth,
    n_reps = nrow(raw),
    fdr_raw = mean(raw$fdp), fdr_rarefied = mean(rar$fdp[i]),
    delta_fdp_mean = mean(d_fdp), delta_fdp_sd = sd(d_fdp),
    delta_fdp_ci_lo = ci_f[1], delta_fdp_ci_hi = ci_f[2],
    delta_fdp_excludes_zero = ci_f[1] > 0 | ci_f[2] < 0,
    frac_reps_delta_pos = mean(d_fdp > 0), frac_reps_delta_zero = mean(d_fdp == 0),
    power_raw = mean(raw$power_itt), power_rarefied = mean(rar$power_itt[i]),
    delta_power_mean = mean(d_pow),
    delta_power_ci_lo = ci_p[1], delta_power_ci_hi = ci_p[2],
    stringsAsFactors = FALSE)
}
res <- do.call(rbind, out)
rownames(res) <- NULL
res <- res[order(res$method, res$cell_id), ]
write.table(res, file.path(out_dir, "b4_delta_fdp_paired_ci.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("done: ", nrow(res), " cell x method rows; bootstrap B = ", BOOT)
print(res[, c("cell_id", "method", "delta_fdp_mean", "delta_fdp_ci_lo",
              "delta_fdp_ci_hi", "delta_fdp_excludes_zero")])
