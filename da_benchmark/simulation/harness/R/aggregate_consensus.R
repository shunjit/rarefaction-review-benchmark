#!/usr/bin/env Rscript
# Consensus (multi-method agreement) metrics from a D-sim run — R3-4.
#
# Reviewer #3 (comment 4, manuscript L233): lowering false positives by
# requiring cross-method agreement is clear, but how can multi-method
# approaches improve false-NEGATIVE rates? This aggregation quantifies both
# directions from the same run:
#   - intersection/majority rules cut FP (FP are method-specific)
#   - the union rule cuts false negatives (methods miss different spiked taxa)
# and records the overlap structure that is the mechanism.
#
# Rules per (cell, rep, condition, taxon), over the three primary methods
# (ancombc2, maaslin3, wilcoxon; maaslin3_prev excluded; boots = its own
# >50% stability call, no extra thresholds here):
#   any1 = union  (>= 1 method detects)
#   maj2 = majority (>= 2 of 3)
#   all3 = intersection (all 3)
# Single-method rows are re-derived under the same pipeline as a consistency
# check against metrics_by_cell.tsv (identical definitions, design v2 s10).
#
# Outputs (run summary/):
#   metrics_consensus_by_cell.tsv  — cell x condition x rule metrics
#   consensus_mechanism_by_cell.tsv — per cell x condition: among union TP
#     (resp. FP), share detected by exactly 1/2/3 methods, plus per-method
#     unique contributions (taxa only that method finds)
# usage: Rscript aggregate_consensus.R <run_dir>

args <- commandArgs(trailingOnly = TRUE)
run_dir <- args[1]
stopifnot(dir.exists(run_dir))
out_dir <- file.path(run_dir, "summary")
dir.create(out_dir, showWarnings = FALSE)

METHODS <- c("ancombc2", "maaslin3", "wilcoxon")
CONDS <- c("raw", "rarefied", "boots")

res_files <- Sys.glob(file.path(run_dir, "results", "*", "rep_[0-9]*.tsv.gz"))
res_files <- res_files[!grepl("_truth", res_files)]
stopifnot(length(res_files) > 0)

# per rep: wide detected matrix (taxon x method) per condition -> long rows
read_one <- function(f) {
  r <- read.delim(gzfile(f), stringsAsFactors = FALSE)
  tf <- sub("\\.tsv\\.gz$", "_truth.tsv.gz", f)
  t <- read.delim(gzfile(tf), stringsAsFactors = FALSE)
  r <- r[r$method %in% METHODS & r$condition %in% CONDS, ]
  r$detected[is.na(r$detected)] <- 0L
  taxa <- t$taxon
  out <- lapply(CONDS, function(cc) {
    rc <- r[r$condition == cc, ]
    det <- sapply(METHODS, function(m) {
      x <- rc[rc$method == m, ]
      x$detected[match(taxa, x$taxon)]
    })
    det[is.na(det)] <- 0L
    data.frame(cell_id = r$cell_id[1], scenario = r$scenario[1],
               rep = r$rep[1], condition = cc, taxon = taxa,
               spiked = t$spiked, tested = t$tested,
               det_a = det[, 1], det_m = det[, 2], det_w = det[, 3],
               n_det = rowSums(det), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}
message(sprintf("reading %d rep files ...", length(res_files)))
d <- do.call(rbind, lapply(res_files, read_one))
d <- d[d$tested == 1L, ]                     # intention-to-test universe
d$is_null <- d$spiked == 0L

# intention-to-test power denominators (spiked_n per cell x rep)
gi_files <- Sys.glob(file.path(run_dir, "results", "*", "rep_[0-9]*_geninfo.tsv"))
gi <- do.call(rbind, lapply(gi_files, read.delim))

RULES <- list(ancombc2 = function(x) x$det_a, maaslin3 = function(x) x$det_m,
              wilcoxon = function(x) x$det_w,
              any1 = function(x) as.integer(x$n_det >= 1L),
              maj2 = function(x) as.integer(x$n_det >= 2L),
              all3 = function(x) as.integer(x$n_det == 3L))

key <- c("cell_id", "scenario", "condition")
grp <- interaction(d[key], drop = TRUE)
summ <- do.call(rbind, lapply(split(d, grp), function(g) {
  do.call(rbind, lapply(names(RULES), function(rn) {
    g$detected <- RULES[[rn]](g)
    per_rep <- split(g, g$rep)
    fdp <- vapply(per_rep, function(r) {
      sum(r$detected & r$is_null) / max(sum(r$detected), 1)
    }, numeric(1))
    fp_n <- vapply(per_rep, function(r) sum(r$detected & r$is_null), numeric(1))
    any_fp <- vapply(per_rep, function(r) any(r$detected & r$is_null), logical(1))
    row <- data.frame(g[1, key], rule = rn,
      n_reps = length(per_rep),
      tested_n_mean = nrow(g) / length(per_rep),
      fpr_taxon = sum(g$detected & g$is_null) / sum(g$is_null),
      mean_fp = mean(fp_n), fwer = mean(any_fp),
      mean_fdp = mean(fdp), fdr_mc_se = sd(fdp) / sqrt(length(per_rep)),
      power_itt = NA_real_, power_cot = NA_real_,
      stringsAsFactors = FALSE)
    if (g$scenario[1] == "spike") {
      tp <- vapply(per_rep, function(r) sum(r$detected & !r$is_null), numeric(1))
      st <- vapply(per_rep, function(r) sum(!r$is_null), numeric(1))
      denom <- gi$spiked_n[match(paste(g$cell_id[1], names(per_rep)),
                                 paste(gi$cell_id, gi$rep))]
      row$power_cot <- mean(tp / pmax(st, 1))
      row$power_itt <- mean(tp / denom)
    }
    row
  }))
}))
rownames(summ) <- NULL
write.table(summ[order(summ$cell_id, summ$condition, summ$rule), ],
            file.path(out_dir, "metrics_consensus_by_cell.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- mechanism: overlap structure of union detections ------------------------
mech <- do.call(rbind, lapply(split(d, grp), function(g) {
  u <- g[g$n_det >= 1L, ]
  tp <- u[!u$is_null, ]; fp <- u[u$is_null, ]
  n_reps <- length(unique(g$rep))
  share <- function(x, k) if (nrow(x)) mean(x$n_det == k) else NA_real_
  uniq <- function(x, col) if (nrow(x)) sum(x[[col]] == 1L & x$n_det == 1L) / n_reps else 0
  data.frame(g[1, key], n_reps = n_reps,
    tp_union_per_rep = nrow(tp) / n_reps, fp_union_per_rep = nrow(fp) / n_reps,
    tp_share_1of3 = share(tp, 1), tp_share_2of3 = share(tp, 2),
    tp_share_3of3 = share(tp, 3),
    fp_share_1of3 = share(fp, 1), fp_share_2of3 = share(fp, 2),
    fp_share_3of3 = share(fp, 3),
    tp_uniq_ancombc2 = uniq(tp, "det_a"), tp_uniq_maaslin3 = uniq(tp, "det_m"),
    tp_uniq_wilcoxon = uniq(tp, "det_w"),
    fp_uniq_ancombc2 = uniq(fp, "det_a"), fp_uniq_maaslin3 = uniq(fp, "det_m"),
    fp_uniq_wilcoxon = uniq(fp, "det_w"),
    stringsAsFactors = FALSE)
}))
rownames(mech) <- NULL
write.table(mech[order(mech$cell_id, mech$condition), ],
            file.path(out_dir, "consensus_mechanism_by_cell.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message(sprintf("consensus metrics over %d cells x %d conditions done",
                length(unique(d$cell_id)), length(CONDS)))
