#!/usr/bin/env Rscript
# D-real per-rep artificial-group construction (D_HARNESS_DESIGN.md v2 §5).
#
# One invocation = one (cell, rep). From the full real count table it draws an
# 80% subsample within each stratum, splits each stratum at its median depth
# (High_c / Low_c), and pairs that depth-ordered assignment with a restricted
# permutation of the SAME subsample (labels shuffled within stratum, group
# sizes fixed). Both assignments share one set of count tables; only the
# metadata differs.
#
# Writes into $WORKDIR:
#   raw.tsv, rarefied.tsv, boots/b###.tsv   (tested taxa x subsample; shared
#       by both assignments — prevalence >= 10% on the raw subsample table)
#   meta_depth.tsv, meta_perm.tsv   sample_id, group, depth_raw, stratum
#       (G1 = Low/shallow, G2 = High/deep — matches the D-sim sign convention)
#   truth.tsv    per-taxon rows over the union universe (count > 0 in the
#       subsample): safe taxon id, original asv_id, empirical pi, abundance
#       tercile stratum, spiked = 0 (no ground truth injected), tested flag
#   geninfo.tsv  one row of rep-level bookkeeping + §5 health checks
#
# Taxon ids are renamed to syntactic ids a00001.. because hex ASV hashes can
# start with a digit, which data.frame/maaslin3 would silently rewrite to
# "X<hash>" and break result matching. The asv_id mapping lives in truth.tsv.
#
# Seeds (§12): subsample+split = dataset_seed; single rarefaction = +901;
# boots replicate b = +b; label permutation = +902.
#
# A-1 (CG-04, 2026-08-16 additive): STRAT_MODE "stratified_extreme" keeps the
# same 80% subsample stream, then keeps only the bottom / top FILTER_VAL
# (e.g. 0.25 = quartile) of each stratum by depth: shallowest -> G1, deepest
# -> G2, middle discarded. Existing modes are byte-identical to before.

source(file.path(Sys.getenv("HARNESS_R_DIR", "."), "lib_common.R"))

REVISION_ROOT <- env_require("REVISION_ROOT")
WORKDIR    <- env_require("WORKDIR")
CELL_ID    <- env_require("CELL_ID")
DATASET    <- env_require("DATASET")
STRAT_COL  <- env_require("STRAT_COL")
STRAT_MODE <- env_require("STRAT_MODE")          # stratified | single | stratified_extreme
FILTER_VAL <- env_require("FILTER_VAL")          # single: stratum kept; stratified_extreme: quantile; NA otherwise
SUB_FRAC   <- env_num("SUBSAMPLE_FRAC")
B          <- as.integer(env_num("B"))
REP        <- as.integer(env_num("REP"))
BASE_SEED  <- env_num("BASE_SEED")
CELL_INDEX <- as.integer(env_num("CELL_INDEX"))

PREV_MIN <- 0.10

stopifnot(STRAT_MODE %in% c("stratified", "single", "stratified_extreme"),
          SUB_FRAC > 0, SUB_FRAC <= 1)
EXT_Q <- NA_real_
if (STRAT_MODE == "stratified_extreme") {
  EXT_Q <- as.numeric(FILTER_VAL)
  stopifnot(is.finite(EXT_Q), EXT_Q > 0, EXT_Q < 0.5)
}

input_dir <- file.path(REVISION_ROOT, "realdata_typeI", "input")
counts_rds <- file.path(input_dir, paste0(DATASET, "_asv_counts.rds"))
meta_path <- file.path(input_dir, paste0(DATASET, "_meta.tsv"))
for (p in c(counts_rds, meta_path)) if (!file.exists(p)) {
  cat(sprintf("ERROR: required input MISSING: %s - stopping.\n", p), file = stderr()); quit(status = 1)
}

meta_all <- read.delim(meta_path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("sample-id", STRAT_COL, "depth") %in% colnames(meta_all)))
meta_all$stratum <- as.character(meta_all[[STRAT_COL]])

pool <- if (STRAT_MODE == "single") meta_all[meta_all$stratum == FILTER_VAL, ] else meta_all
if (nrow(pool) < 8L) {
  cat(sprintf("ERROR: pool has only %d samples - stopping.\n", nrow(pool)), file = stderr()); quit(status = 1)
}

X_full <- readRDS(counts_rds)                    # integer matrix taxa x samples
stopifnot(all(pool[["sample-id"]] %in% colnames(X_full)))

ds <- dataset_seed(BASE_SEED, CELL_INDEX, REP)

# --- 1. 80% subsample within stratum, median-depth split (High_c / Low_c) ----
set.seed(ds)
strata <- sort(unique(pool$stratum))
sub_ids <- character(0); grp_depth <- character(0)
for (st in strata) {
  ids <- pool[["sample-id"]][pool$stratum == st]
  n_keep <- round(SUB_FRAC * length(ids))
  keep <- if (n_keep < length(ids)) sample(ids, n_keep) else ids
  d <- pool$depth[match(keep, pool[["sample-id"]])]
  keep <- keep[order(d, keep)]                   # deterministic tie-break by id
  if (STRAT_MODE == "stratified_extreme") {
    n_q <- floor(EXT_Q * n_keep)                 # per-tail size in this stratum
    if (n_q < 2L) {
      cat(sprintf("ERROR: stratum %s too small for extreme split (n_keep=%d, q=%.2f) - stopping.\n",
                  st, n_keep, EXT_Q), file = stderr()); quit(status = 1)
    }
    keep <- c(keep[seq_len(n_q)], keep[seq(n_keep - n_q + 1L, n_keep)])
    g <- rep(c("G1", "G2"), c(n_q, n_q))         # shallowest tail / deepest tail
  } else {
    n_low <- floor(n_keep / 2)                   # odd n: extra sample goes High
    g <- rep(c("G1", "G2"), c(n_low, n_keep - n_low))
  }
  sub_ids <- c(sub_ids, keep); grp_depth <- c(grp_depth, g)
}
names(grp_depth) <- sub_ids

# --- 2. paired control: restricted permutation within stratum ---------------
set.seed(ds + SEED_OFF_PERMUTE)
grp_perm <- grp_depth
strat_of <- pool$stratum[match(sub_ids, pool[["sample-id"]])]
names(strat_of) <- sub_ids
for (st in strata) {
  ids <- sub_ids[strat_of == st]
  grp_perm[ids] <- sample(grp_depth[ids])        # sizes fixed by construction
}

# --- 3. subset counts, safe taxon ids, tested set ----------------------------
X <- X_full[, sub_ids, drop = FALSE]
rm(X_full)
X <- X[rowSums(X) > 0L, , drop = FALSE]          # union universe of the subsample
asv_ids <- rownames(X)
taxa <- sprintf("a%05d", seq_len(nrow(X)))
rownames(X) <- taxa

prevalence_raw <- rowMeans(X > 0)
tested <- prevalence_raw >= PREV_MIN
depth_raw <- colSums(X)

# --- 4. rarefied + boots tables (sparse-aware; zeros carry no information) ---
rarefy_vec_sparse <- function(x, depth) {
  idx <- which(x > 0L)
  out <- integer(length(x))
  out[idx] <- rarefy_vec(x[idx], depth)
  out
}
rarefy_mat_sparse <- function(m, depth) {
  r <- vapply(seq_len(ncol(m)), function(j) rarefy_vec_sparse(m[, j], depth), integer(nrow(m)))
  dimnames(r) <- dimnames(m)
  r
}

raref_depth <- min(depth_raw)                    # zero-dropout default (§2)
set.seed(ds + SEED_OFF_SINGLE_RAREF)
R1 <- rarefy_mat_sparse(X, raref_depth)
dir.create(file.path(WORKDIR, "boots"), recursive = TRUE, showWarnings = FALSE)
for (b in seq_len(B)) {
  set.seed(ds + b)
  write_counts(rarefy_mat_sparse(X, raref_depth)[tested, , drop = FALSE],
               file.path(WORKDIR, "boots", sprintf("b%03d.tsv", b)))
}
write_counts(X[tested, , drop = FALSE], file.path(WORKDIR, "raw.tsv"))
write_counts(R1[tested, , drop = FALSE], file.path(WORKDIR, "rarefied.tsv"))

for (v in c("depth", "perm")) {
  g <- if (v == "depth") grp_depth else grp_perm
  write.table(data.frame(sample_id = sub_ids, group = g[sub_ids],
                         depth_raw = depth_raw[sub_ids], stratum = strat_of[sub_ids]),
              file.path(WORKDIR, sprintf("meta_%s.tsv", v)),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# --- 5. truth over the union universe (no injected signal: spiked = 0) -------
pi_emp <- rowMeans(sweep(X, 2, pmax(depth_raw, 1L), "/"))
strat_tax <- rep(c("high", "mid", "low"),
                 times = diff(round(seq(0, nrow(X), length.out = 4))))[order(order(pi_emp, decreasing = TRUE))]
write.table(data.frame(taxon = taxa, asv_id = asv_ids, pi = pi_emp,
                       stratum = strat_tax, spiked = 0L, direction = NA_character_,
                       mult_absolute = 1, rel_fc_analytic = 1,
                       tested = as.integer(tested), prevalence_raw = prevalence_raw),
            file.path(WORKDIR, "truth.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# --- 6. geninfo: bookkeeping + §5 health checks -------------------------------
med <- function(ids) median(depth_raw[ids])
g1_d <- sub_ids[grp_depth[sub_ids] == "G1"]; g2_d <- sub_ids[grp_depth[sub_ids] == "G2"]
g1_p <- sub_ids[grp_perm[sub_ids] == "G1"]; g2_p <- sub_ids[grp_perm[sub_ids] == "G2"]
comp <- vapply(strata, function(st) {
  ids <- sub_ids[strat_of == st]
  sprintf("%s:%d/%d", st, sum(grp_depth[ids] == "G1"), sum(grp_depth[ids] == "G2"))
}, character(1))
comp_perm_equal <- all(vapply(strata, function(st) {
  ids <- sub_ids[strat_of == st]
  sum(grp_perm[ids] == "G1") == sum(grp_depth[ids] == "G1")
}, logical(1)))

write.table(data.frame(cell_id = CELL_ID, scenario = "real", rep = REP,
                       dataset = DATASET, strat_col = STRAT_COL, strat_mode = STRAT_MODE,
                       filter_val = FILTER_VAL, subsample_frac = SUB_FRAC,
                       n_pool = nrow(pool), n_sub = length(sub_ids),
                       n_g1 = length(g1_d), n_g2 = length(g2_d),
                       composition_g1_g2 = paste(comp, collapse = "|"),
                       comp_perm_equal = comp_perm_equal,
                       median_depth_g1 = med(g1_d), median_depth_g2 = med(g2_d),
                       achieved_ratio_depth = med(g2_d) / med(g1_d),
                       achieved_ratio_perm = med(g2_p) / med(g1_p),
                       raref_depth = raref_depth, raref_dropout_n = 0L,
                       union_taxa_n = nrow(X), tested_n = sum(tested),
                       zero_frac_raw = mean(X[tested, ] == 0),
                       zero_frac_rarefied = mean(R1[tested, ] == 0),
                       B = B, base_seed = BASE_SEED, cell_index = CELL_INDEX,
                       dataset_seed = ds),
            file.path(WORKDIR, "geninfo.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("generated %s rep %d: %d samples (G1 %d / G2 %d), tested %d/%d taxa, raref depth %d, achieved rho %.2f (perm %.2f)\n",
            CELL_ID, REP, length(sub_ids), length(g1_d), length(g2_d),
            sum(tested), nrow(X), raref_depth,
            med(g2_d) / med(g1_d), med(g2_p) / med(g1_p)))
