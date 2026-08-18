#!/usr/bin/env Rscript
# B-1 (CG-05, GM-01): ANCOM-BC2 diff_robust re-aggregation. No refit —
# pure re-aggregation of the frozen per-rep mirrors.
#   default     = detected               (primary rule, q_group < 0.05)
#   diff_robust = detected & passed_ss   (also survives the pseudo-count
#                 sensitivity filter; saved for raw/rarefied only — the
#                 pseudo_sens hybrid of 判断2 leaves boots out by design)
# The default-rule output is asserted against the run's frozen summary
# (summary/metrics_by_cell.tsv) as a self-check before anything is written.
# usage: Rscript b1_diff_robust.R <dsim_run_dir> <dreal_run_dir> <out_dir>

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
dsim_dir  <- args[1]; dreal_dir <- args[2]; out_dir <- args[3]
stopifnot(dir.exists(dsim_dir), dir.exists(dreal_dir))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

list_rep_files <- function(run_dir) {
  f <- Sys.glob(file.path(run_dir, "results", "*", "rep_[0-9]*.tsv.gz"))
  f <- f[!grepl("_truth", f)]
  stopifnot(length(f) > 0)
  f
}

read_run <- function(run_dir) {
  files <- list_rep_files(run_dir)
  message(sprintf("reading %d rep files from %s ...", length(files), basename(run_dir)))
  lst <- lapply(files, function(f) {
    r <- read.delim(gzfile(f), stringsAsFactors = FALSE)
    r <- r[r$method == "ancombc2" & r$condition %in% c("raw", "rarefied"),
           c("cell_id", "scenario", "rep", "condition", "taxon",
             "detected", "reason", "passed_ss")]
    tf <- sub("\\.tsv\\.gz$", "_truth.tsv.gz", f)
    t <- read.delim(gzfile(tf), stringsAsFactors = FALSE)
    if (!"spiked" %in% names(t)) t$spiked <- 0L
    merge(r, t[, c("taxon", "spiked", "tested")], by = "taxon", sort = FALSE)
  })
  d <- do.call(rbind, lst)
  d <- d[d$tested == 1L, ]                    # intention-to-test universe
  d$detected[is.na(d$detected)] <- 0L
  d$passed_ss[is.na(d$passed_ss)] <- 0L
  d$det_default <- as.integer(d$detected == 1L)
  d$det_robust  <- as.integer(d$detected == 1L & d$passed_ss == 1L)
  d$is_null <- d$spiked == 0L
  d
}

dsim  <- read_run(dsim_dir)
dreal <- read_run(dreal_dir)

gi_sim <- do.call(rbind, lapply(
  Sys.glob(file.path(dsim_dir, "results", "*", "rep_[0-9]*_geninfo.tsv")), read.delim))
gi_real <- do.call(rbind, lapply(
  Sys.glob(file.path(dreal_dir, "results", "*", "rep_[0-9]*_geninfo.tsv")), read.delim))

# --- D-sim: same metric construction as aggregate_run.R, per rule ------------
dsim_metrics <- function(d, det_col) {
  key <- c("cell_id", "scenario", "condition")
  out <- do.call(rbind, lapply(split(d, interaction(d[key], drop = TRUE)), function(g) {
    per_rep <- split(g, g$rep)
    det <- function(r) r[[det_col]]
    fdp    <- vapply(per_rep, function(r) sum(det(r) & r$is_null) / max(sum(det(r)), 1), numeric(1))
    fp_n   <- vapply(per_rep, function(r) sum(det(r) & r$is_null), numeric(1))
    any_fp <- vapply(per_rep, function(r) any(det(r) & r$is_null), logical(1))
    row <- data.frame(g[1, key],
      n_reps = length(per_rep),
      tested_n_mean = nrow(g) / length(per_rep),
      fpr_taxon = sum(g[[det_col]] & g$is_null) / sum(g$is_null),
      mean_fp = mean(fp_n), fwer = mean(any_fp),
      mean_fdp = mean(fdp), fdr_mc_se = sd(fdp) / sqrt(length(per_rep)),
      power_itt = NA_real_, power_cot = NA_real_,
      stringsAsFactors = FALSE)
    if (g$scenario[1] == "spike") {
      tp <- vapply(per_rep, function(r) sum(det(r) & !r$is_null), numeric(1))
      st <- vapply(per_rep, function(r) sum(!r$is_null), numeric(1))
      row$power_cot <- mean(tp / pmax(st, 1))
      denom <- gi_sim$spiked_n[match(paste(g$cell_id[1], names(per_rep)),
                                     paste(gi_sim$cell_id, gi_sim$rep))]
      row$power_itt <- mean(tp / denom)
    }
    row
  }))
  rownames(out) <- NULL
  out
}

m_def <- dsim_metrics(dsim, "det_default")
m_rob <- dsim_metrics(dsim, "det_robust")

# self-check: default rule must reproduce the frozen summary exactly
ref <- read.delim(file.path(dsim_dir, "summary", "metrics_by_cell.tsv"))
ref <- ref[ref$method == "ancombc2" & ref$condition %in% c("raw", "rarefied"), ]
chk <- merge(ref, m_def, by = c("cell_id", "scenario", "condition"),
             suffixes = c("_ref", "_new"))
stopifnot(nrow(chk) == nrow(m_def))
for (m in c("fpr_taxon", "mean_fp", "fwer", "mean_fdp", "fdr_mc_se",
            "power_itt", "power_cot")) {
  dv <- abs(chk[[paste0(m, "_ref")]] - chk[[paste0(m, "_new")]])
  dmax <- suppressWarnings(max(dv, na.rm = TRUE))
  if (is.finite(dmax)) {
    message(sprintf("self-check dsim default %-9s: max |diff| = %.3g", m, dmax))
    stopifnot(dmax < 1e-9)
  }
}

m_def$rule <- "default"; m_rob$rule <- "diff_robust"
long <- rbind(m_def, m_rob)
write.table(long[order(long$cell_id, long$condition, long$rule), ],
            file.path(out_dir, "b1_dsim_ancombc2_rules_by_cell.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

met <- c("fpr_taxon", "mean_fp", "fwer", "mean_fdp", "fdr_mc_se",
         "power_itt", "power_cot")
wide <- merge(m_def[, c("cell_id", "scenario", "condition", "n_reps",
                        "tested_n_mean", met)],
              m_rob[, c("cell_id", "scenario", "condition", met)],
              by = c("cell_id", "scenario", "condition"),
              suffixes = c("_default", "_robust"))
tot <- aggregate(cbind(det_default, det_robust) ~ cell_id + scenario + condition,
                 data = dsim, FUN = sum)
names(tot)[names(tot) == "det_default"] <- "n_det_total_default"
names(tot)[names(tot) == "det_robust"]  <- "n_det_total_robust"
tot$ss_survival_of_detections <- ifelse(tot$n_det_total_default > 0,
  tot$n_det_total_robust / tot$n_det_total_default, NA_real_)
wide <- merge(wide, tot, by = c("cell_id", "scenario", "condition"))
write.table(wide[order(wide$cell_id, wide$condition), ],
            file.path(out_dir, "b1_dsim_ancombc2_rules_wide.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- D-real: descriptive re-aggregation per rule (aggregate_real_run.R form) --
real_per_rep <- function(d, det_col) {
  key <- c("cell_id", "condition", "scenario")
  out <- do.call(rbind, lapply(split(d, interaction(d[c(key, "rep")], drop = TRUE)),
    function(g) data.frame(g[1, key], rep = g$rep[1],
      n_tested = nrow(g), n_detected = sum(g[[det_col]]),
      any_detected = as.integer(any(g[[det_col]] == 1L)),
      stringsAsFactors = FALSE)))
  rownames(out) <- NULL
  out
}
pr_def <- real_per_rep(dreal, "det_default")
pr_rob <- real_per_rep(dreal, "det_robust")

# self-check against frozen dreal summary (median/any-rate for ancombc2)
q1 <- function(x) quantile(x, 0.25, names = FALSE)
q3 <- function(x) quantile(x, 0.75, names = FALSE)
real_summ <- function(pr) {
  key <- c("cell_id", "condition", "scenario")
  out <- do.call(rbind, lapply(split(pr, interaction(pr[key], drop = TRUE)),
    function(g) data.frame(g[1, key], n_reps = nrow(g),
      tested_n_median = median(g$n_tested),
      det_median = median(g$n_detected), det_q25 = q1(g$n_detected),
      det_q75 = q3(g$n_detected), det_min = min(g$n_detected),
      det_max = max(g$n_detected),
      any_detection_rate = mean(g$any_detected),
      stringsAsFactors = FALSE)))
  rownames(out) <- NULL
  out
}
rs_def <- real_summ(pr_def)
rs_rob <- real_summ(pr_rob)

ref_r <- read.delim(file.path(dreal_dir, "summary", "metrics_real_by_cell.tsv"))
ref_r <- ref_r[ref_r$method == "ancombc2" & ref_r$condition %in% c("raw", "rarefied"), ]
chk_r <- merge(ref_r, rs_def, by = c("cell_id", "condition", "scenario"),
               suffixes = c("_ref", "_new"))
stopifnot(nrow(chk_r) == nrow(rs_def))
for (m in c("det_median", "det_q25", "det_q75", "any_detection_rate")) {
  dmax <- max(abs(chk_r[[paste0(m, "_ref")]] - chk_r[[paste0(m, "_new")]]))
  message(sprintf("self-check dreal default %-18s: max |diff| = %.3g", m, dmax))
  stopifnot(dmax < 1e-9)
}

rs_def$rule <- "default"; rs_rob$rule <- "diff_robust"
rlong <- rbind(rs_def, rs_rob)
write.table(rlong[order(rlong$cell_id, rlong$condition, rlong$scenario, rlong$rule), ],
            file.path(out_dir, "b1_dreal_ancombc2_rules_by_cell.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# combined per-rep table (both rules side by side + achieved depth ratios)
key_r <- c("cell_id", "condition", "scenario", "rep")
prw <- merge(pr_def[, c(key_r, "n_tested", "n_detected", "any_detected")],
             pr_rob[, c(key_r, "n_detected", "any_detected")],
             by = key_r, suffixes = c("_default", "_robust"))
m <- match(paste(prw$cell_id, prw$rep), paste(gi_real$cell_id, gi_real$rep))
prw$achieved_ratio_depth <- gi_real$achieved_ratio_depth[m]
prw$achieved_ratio_perm  <- gi_real$achieved_ratio_perm[m]
prw$raref_depth <- gi_real$raref_depth[m]
prw$n_sub <- gi_real$n_sub[m]
write.table(prw[order(prw$cell_id, prw$condition, prw$scenario, prw$rep), ],
            file.path(out_dir, "b1_dreal_ancombc2_rules_per_rep.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# paired within-rep differences (depth - perm), per rule
paired_summ <- function(pr, rule_name) {
  wk <- c("cell_id", "condition", "rep")
  dep <- pr[pr$scenario == "real_depth", ]
  per <- pr[pr$scenario == "real_perm", ]
  mm <- match(interaction(dep[wk], drop = FALSE), interaction(per[wk], drop = FALSE))
  pair <- dep[!is.na(mm), wk]
  pair$diff_detected <- dep$n_detected[!is.na(mm)] - per$n_detected[mm[!is.na(mm)]]
  out <- do.call(rbind, lapply(
    split(pair, interaction(pair[c("cell_id", "condition")], drop = TRUE)),
    function(g) data.frame(g[1, c("cell_id", "condition")], n_pairs = nrow(g),
      diff_median = median(g$diff_detected), diff_q25 = q1(g$diff_detected),
      diff_q75 = q3(g$diff_detected),
      frac_pairs_depth_higher = mean(g$diff_detected > 0),
      frac_pairs_equal = mean(g$diff_detected == 0),
      stringsAsFactors = FALSE)))
  rownames(out) <- NULL
  out$rule <- rule_name
  out
}
pair_both <- rbind(paired_summ(pr_def, "default"), paired_summ(pr_rob, "diff_robust"))
write.table(pair_both[order(pair_both$cell_id, pair_both$condition, pair_both$rule), ],
            file.path(out_dir, "b1_dreal_ancombc2_rules_paired.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- provenance ---------------------------------------------------------------
run_id_of <- function(run_dir) {
  y <- readLines(file.path(run_dir, "config.yaml"), warn = FALSE)
  sub("^run_id:\\s*", "", y[grepl("^run_id:", y)][1])
}
prov <- c(
  "B-1 diff_robust re-aggregation (CG-05, GM-01) — no refit",
  sprintf("generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("script: b1_diff_robust.R"),
  sprintf("dsim run: %s (%d rep files, %d cells)", run_id_of(dsim_dir),
          length(list_rep_files(dsim_dir)), length(unique(dsim$cell_id))),
  sprintf("dreal run: %s (%d rep files, %d cells)", run_id_of(dreal_dir),
          length(list_rep_files(dreal_dir)), length(unique(dreal$cell_id))),
  "rules: default = detected (q<0.05); diff_robust = detected & passed_ss",
  "scope: method = ancombc2, condition in {raw, rarefied} (boots has pseudo_sens=FALSE by design)",
  "self-check: default-rule metrics reproduced the frozen run summaries to < 1e-9"
)
writeLines(prov, file.path(out_dir, "b1_provenance.txt"))

message("done. outputs in ", out_dir)
print(wide[wide$cell_id %in% c("null_r1", "null_r10", "null_r100",
                               "spike_r1_f3", "spike_r10_f3", "spike_r100_f3"),
           c("cell_id", "condition", "fwer_default", "fwer_robust",
             "mean_fdp_default", "mean_fdp_robust", "power_itt_default",
             "power_itt_robust", "ss_survival_of_detections")])
