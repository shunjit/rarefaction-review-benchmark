#!/usr/bin/env Rscript
# D-sim calibration step 2: Rice-calibrated Dirichlet-multinomial parameters.
#
# Input : calibrate/output/rice_rhizoplane_asv_counts.tsv (features x samples)
# Output: pi_hat  = empirical mean relative abundance over the top-K ASVs
#                   (renormalised per sample; DM-MLE alpha kept as reference —
#                   the MLE shrinks rare taxa upward and understates sparsity)
#         theta_hat = calibrated by matching the simulated zero fraction to the
#                   real one at a common depth DC (bisection); the DM MLE theta
#                   from dmn(k = 1) is recorded alongside as reference.
#                   Zero fraction is the tuned statistic because rare-taxon
#                   dropout is the phenomenon under study; variance/mean,
#                   Shannon and prevalence remain untuned reality checks.
#
# Design reference: D_HARNESS_DESIGN.md v2 section 4.
# Run: ~/miniforge3/envs/qiime2-amplicon-2025.10/bin/Rscript fit_dm.R <outdir>

suppressMessages(library(DirichletMultinomial))

K <- 300L
DC <- 5000L          # common depth for calibration target + reality check
FIT_SEED <- 20260811
GEN_SEED <- 20260812
ZF_REPS <- 3L        # zero-fraction evaluations averaged per theta candidate

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
stopifnot(!is.na(outdir), dir.exists(outdir))
infile <- file.path(outdir, "rice_rhizoplane_asv_counts.tsv")
stopifnot(file.exists(infile))

counts <- as.matrix(read.delim(infile, row.names = 1, check.names = FALSE))
storage.mode(counts) <- "integer"
cat(sprintf("input: %d ASVs x %d samples\n", nrow(counts), ncol(counts)))

feature_totals <- rowSums(counts)
ord <- order(feature_totals, decreasing = TRUE)
top_idx <- ord[seq_len(K)]
top <- counts[top_idx, , drop = FALSE]
tail_col <- colSums(counts[-top_idx, , drop = FALSE])
coverage <- sum(top) / sum(counts)

# --- pi: empirical mean relative abundance within the top-K subtable ---------
top_sums <- colSums(top)
stopifnot(all(top_sums > 0))
rel <- sweep(top, 2, top_sums, "/")
pi_hat <- rowMeans(rel)
pi_hat <- pi_hat / sum(pi_hat)

# --- reference DM MLE (theta and alpha, for the record) ----------------------
mat <- t(rbind(top, TAIL = tail_col))
fit <- dmn(mat, k = 1, verbose = FALSE, seed = FIT_SEED)
theta_mle <- mixturewt(fit)$theta
alpha <- fitted(fit)[, 1]
pi_mle <- alpha[seq_len(K)] / sum(alpha[seq_len(K)])

# --- real comparator: top-K subtable rarefied to DC --------------------------
rarefy_vec <- function(x, depth) {
  pool <- rep.int(seq_along(x), x)
  tabulate(sample(pool, depth), nbins = length(x))
}
shannon <- function(x) { p <- x[x > 0] / sum(x); -sum(p * log(p)) }

set.seed(GEN_SEED)
keep <- which(top_sums >= DC)
real_rare <- vapply(keep, function(j) rarefy_vec(top[, j], DC), integer(K))
n_gen <- length(keep)
real_zf <- mean(real_rare == 0)

# --- calibrate theta by matching the zero fraction (bisection, monotone) -----
sim_matrix <- function(theta, n, seed) {
  set.seed(seed)
  g <- matrix(rgamma(n * K, shape = theta * pi_hat), nrow = K)
  p <- sweep(g, 2, colSums(g), "/")
  vapply(seq_len(n), function(j) rmultinom(1, DC, p[, j])[, 1], integer(K))
}
zf_at <- function(theta) {
  mean(vapply(seq_len(ZF_REPS),
              function(r) mean(sim_matrix(theta, n_gen, GEN_SEED + 100L * r) == 0),
              numeric(1)))
}
lo <- 5; hi <- 500
stopifnot(zf_at(lo) > real_zf, zf_at(hi) < real_zf)
for (i in seq_len(40)) {
  mid <- sqrt(lo * hi)
  if (zf_at(mid) > real_zf) lo <- mid else hi <- mid
  if (hi / lo < 1.005) break
}
theta_cal <- sqrt(lo * hi)
achieved_zf <- zf_at(theta_cal)

write.table(data.frame(asv_id = rownames(top), pi = pi_hat, pi_dm_mle = pi_mle),
            file.path(outdir, "pi_hat.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(theta_hat = theta_cal, theta_dm_mle = theta_mle,
                       calibration_target = "zero_fraction", target_zf = real_zf,
                       achieved_zf = achieved_zf, depth_dc = DC, K = K,
                       n_samples = ncol(counts), n_samples_zf = n_gen,
                       read_coverage_topK = coverage,
                       fit_seed = FIT_SEED, gen_seed = GEN_SEED),
            file.path(outdir, "theta_hat.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("theta_cal = %.2f (DM MLE reference %.2f), target zf %.3f -> achieved %.3f\n",
            theta_cal, theta_mle, real_zf, achieved_zf))

# --- reality check at (pi_hat, theta_cal) ------------------------------------
sim <- sim_matrix(theta_cal, n_gen, GEN_SEED + 9999L)
per_taxon <- function(m, label) data.frame(
  taxon_rank = seq_len(K), source = label,
  prevalence = rowMeans(m > 0), mean_count = rowMeans(m), var_count = apply(m, 1, var),
  mean_rel = rowMeans(sweep(m, 2, colSums(m), "/")))
pt <- rbind(per_taxon(real_rare, "real"), per_taxon(sim, "sim"))
write.table(pt, file.path(outdir, "reality_per_taxon.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
ps <- data.frame(sample_index = rep(seq_len(n_gen), 2),
                 source = rep(c("real", "sim"), each = n_gen),
                 shannon = c(apply(real_rare, 2, shannon), apply(sim, 2, shannon)))
write.table(ps, file.path(outdir, "reality_per_sample.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

qs <- function(x) sprintf("%.3f / %.3f / %.3f", quantile(x, .25), median(x), quantile(x, .75))
summary_md <- c(
  "# Reality check — Rice-calibrated DM generator vs Rhizoplane top-300 (both at depth 5,000)",
  "",
  sprintf("- parameters: pi = empirical mean composition (top-%d, coverage %.3f of reads); theta = %.1f, calibrated to match the real zero fraction (DM MLE reference: %.1f — shrinks rare taxa and understates sparsity, kept for the record)",
          K, coverage, theta_cal, theta_mle),
  sprintf("- samples compared: %d (real samples with >= %d reads on top-%d ASVs)", n_gen, DC, K),
  sprintf("- zero fraction: real %.3f vs sim %.3f  [tuned statistic]", real_zf, mean(sim == 0)),
  sprintf("- Shannon (Q1/median/Q3): real %s vs sim %s  [untuned]",
          qs(ps$shannon[ps$source == "real"]), qs(ps$shannon[ps$source == "sim"])),
  sprintf("- prevalence, rank-matched real-vs-sim Spearman: %.3f  [untuned]",
          cor(pt$prevalence[pt$source == "real"], pt$prevalence[pt$source == "sim"], method = "spearman")),
  sprintf("- log10 mean count, real-vs-sim Pearson: %.3f  [near-construction: pi is the empirical mean]",
          cor(log10(pt$mean_count[pt$source == "real"] + 0.5), log10(pt$mean_count[pt$source == "sim"] + 0.5))),
  sprintf("- mean-variance: real median var/mean %.1f vs sim %.1f (taxa with mean >= 1)  [untuned]",
          median((pt$var_count / pt$mean_count)[pt$source == "real" & pt$mean_count >= 1]),
          median((pt$var_count / pt$mean_count)[pt$source == "sim" & pt$mean_count >= 1])),
  "",
  "Known limitations (recorded for the response letter): a single-component DM cannot",
  "reproduce the full between-sample heterogeneity of Rhizoplane (field vs greenhouse",
  "runs), so the simulated Shannon spread is narrower than the real one, taxon-level",
  "prevalence ranks decouple from mean-abundance ranks in real data but not under DM",
  "(hence the modest prevalence Spearman), and matching the zero fraction overshoots",
  "the real var/mean ratio by ~1.5x. Mean abundance and global sparsity were",
  "prioritised because effect sizes and rarefaction dropout act on those axes; the",
  "model family (Dirichlet + multinomial) is the one specified by Reviewer #2.",
  "",
  "Per-taxon and per-sample tables: reality_per_taxon.tsv / reality_per_sample.tsv.",
  sprintf("Generator: p ~ Dirichlet(theta * pi), x ~ Multinomial(%d, p); seeds %d/%d.",
          DC, FIT_SEED, GEN_SEED))
writeLines(summary_md, file.path(outdir, "reality_check.md"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_fit_dm.txt"))
cat(paste(summary_md[3:9], collapse = "\n"), "\n")
