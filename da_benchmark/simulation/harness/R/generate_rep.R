#!/usr/bin/env Rscript
# D-sim per-rep data generation (D_HARNESS_DESIGN.md v2 sections 4, 6, 10, 12).
#
# One invocation = one (cell, rep). Writes into $WORKDIR:
#   raw.tsv, rarefied.tsv, boots/b###.tsv  (taxa x samples, subset to the
#       intention-to-test set = prevalence >= 10% on the raw table)
#   meta.tsv     sample_id, group, depth_raw
#   truth.tsv    per-taxon ground truth incl. tested flag (ALL K taxa)
#   geninfo.tsv  one row of rep-level bookkeeping
#
# Spike-in happens on the LATENT scale: p ~ Dirichlet(theta*pi) is treated as
# latent absolute abundance, group-2 spiked taxa are multiplied by f or 1/f,
# and multinomial sampling renormalises. Non-spiked taxa are absolute nulls.

source(file.path(Sys.getenv("HARNESS_R_DIR", "."), "lib_common.R"))

REVISION_ROOT <- env_require("REVISION_ROOT")
WORKDIR   <- env_require("WORKDIR")
CELL_ID   <- env_require("CELL_ID")
SCENARIO  <- env_require("SCENARIO")            # null | spike
K         <- as.integer(env_num("K"))
N_GRP     <- as.integer(env_num("N_PER_GROUP"))
RHO       <- env_num("RHO")
D0        <- env_num("D0")
THETA_RAW <- env_require("THETA")               # numeric or "cal"
FC        <- env_num("FC")                      # 1 when null
DA_FRAC   <- env_num("DA_FRAC")                 # 0 when null
B         <- as.integer(env_num("B"))
REP       <- as.integer(env_num("REP"))
BASE_SEED <- env_num("BASE_SEED")
CELL_INDEX<- as.integer(env_num("CELL_INDEX"))

PREV_MIN  <- 0.10
DEPTH_FLOOR <- 500

calib_dir <- file.path(REVISION_ROOT, "simulation", "harness", "calibrate", "output")
pi_tab <- read.delim(file.path(calib_dir, "pi_hat.tsv"))
stopifnot(nrow(pi_tab) == K)
pi_hat <- pi_tab$pi
theta <- if (THETA_RAW == "cal") read.delim(file.path(calib_dir, "theta_hat.tsv"))$theta_hat else as.numeric(THETA_RAW)
stopifnot(is.finite(theta), theta > 0)

ds <- dataset_seed(BASE_SEED, CELL_INDEX, REP)
set.seed(ds)

n_tot <- 2L * N_GRP
group <- rep(c("G1", "G2"), each = N_GRP)
taxa <- sprintf("t%03d", seq_len(K))
samples <- sprintf("s%03d", seq_len(n_tot))

# --- 1. depths (G2 median = D0 * rho) ---------------------------------------
depths <- pmax(DEPTH_FLOOR, round(exp(rnorm(n_tot, log(D0), 0.3))))
depths[group == "G2"] <- pmax(DEPTH_FLOOR, round(exp(rnorm(N_GRP, log(D0 * RHO), 0.3))))

# --- 2. spike multipliers (latent scale, stratified by abundance terciles) ---
mult <- rep(1, K); spiked <- rep(0L, K); direction <- rep(NA_character_, K)
if (SCENARIO == "spike") {
  n_spike <- round(DA_FRAC * K)
  strat <- rep(c("high", "mid", "low"), each = ceiling(K / 3))[order(order(pi_hat, decreasing = TRUE))]
  per_str <- diff(round(seq(0, n_spike, length.out = 4)))    # near-equal split over 3 strata
  chosen <- unlist(lapply(seq_along(c("high", "mid", "low")), function(si) {
    pool <- which(strat == c("high", "mid", "low")[si])
    sample(pool, per_str[si])
  }))
  dirs <- sample(rep(c("up", "down"), length.out = length(chosen)))
  mult[chosen] <- ifelse(dirs == "up", FC, 1 / FC)
  spiked[chosen] <- 1L; direction[chosen] <- dirs
}
strat_all <- rep(c("high", "mid", "low"), each = ceiling(K / 3))[order(order(pi_hat, decreasing = TRUE))]
rel_fc_analytic <- mult / sum(pi_hat * mult)     # first-order realized relative FC

# --- 3. latent composition + counts ------------------------------------------
X <- matrix(0L, nrow = K, ncol = n_tot, dimnames = list(taxa, samples))
for (i in seq_len(n_tot)) {
  a <- rgamma(K, shape = theta * pi_hat)
  a[a == 0] <- .Machine$double.xmin              # guard against exact-0 gamma draws
  if (group[i] == "G2") a <- a * mult
  X[, i] <- rmultinom(1, depths[i], a / sum(a))[, 1]
}

# --- 4. intention-to-test set (raw prevalence >= 10%, common to all arms) ----
prevalence_raw <- rowMeans(X > 0)
tested <- prevalence_raw >= PREV_MIN

# --- 5. rarefied + boots tables (full matrix first, then subset) -------------
raref_depth <- min(colSums(X))
set.seed(ds + SEED_OFF_SINGLE_RAREF)
R1 <- rarefy_mat(X, raref_depth)
dir.create(file.path(WORKDIR, "boots"), recursive = TRUE, showWarnings = FALSE)
for (b in seq_len(B)) {
  set.seed(ds + b)
  write_counts(rarefy_mat(X, raref_depth)[tested, , drop = FALSE],
               file.path(WORKDIR, "boots", sprintf("b%03d.tsv", b)))
}
write_counts(X[tested, , drop = FALSE], file.path(WORKDIR, "raw.tsv"))
write_counts(R1[tested, , drop = FALSE], file.path(WORKDIR, "rarefied.tsv"))

write.table(data.frame(sample_id = samples, group = group, depth_raw = colSums(X)),
            file.path(WORKDIR, "meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(taxon = taxa, pi = pi_hat, stratum = strat_all,
                       spiked = spiked, direction = direction, mult_absolute = mult,
                       rel_fc_analytic = rel_fc_analytic, tested = as.integer(tested),
                       prevalence_raw = prevalence_raw),
            file.path(WORKDIR, "truth.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(cell_id = CELL_ID, scenario = SCENARIO, rep = REP,
                       K = K, n_per_group = N_GRP, rho = RHO, d0 = D0,
                       theta = theta, theta_mode = THETA_RAW, fc = FC, da_frac = DA_FRAC,
                       B = B, base_seed = BASE_SEED, cell_index = CELL_INDEX,
                       dataset_seed = ds,
                       median_depth_g1 = median(depths[group == "G1"]),
                       median_depth_g2 = median(depths[group == "G2"]),
                       achieved_rho = median(depths[group == "G2"]) / median(depths[group == "G1"]),
                       raref_depth = raref_depth, tested_n = sum(tested),
                       zero_frac_raw = mean(X == 0), zero_frac_rarefied = mean(R1 == 0),
                       spiked_n = sum(spiked), spiked_tested_n = sum(spiked == 1L & tested)),
            file.path(WORKDIR, "geninfo.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("generated %s rep %d: tested %d/%d taxa, raref depth %d, achieved rho %.1f\n",
            CELL_ID, REP, sum(tested), K, raref_depth,
            median(depths[group == "G2"]) / median(depths[group == "G1"])))
