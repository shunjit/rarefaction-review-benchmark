#!/usr/bin/env Rscript
# ANCOM-BC2 runner (D_HARNESS_DESIGN.md v2 section 11).
# Fixed config: fix_formula="group", p_adj_method="BH", prv_cut=0 (global
# pre-filter already applied), lib_cut=0, struc_zero=FALSE, neg_lb=FALSE,
# alpha=0.05, global=FALSE, pairwise=FALSE, n_cl=1.
# Primary call: q_groupG2 < 0.05. passed_ss recorded as secondary column.
# PSEUDO_SENS env (TRUE/FALSE, default TRUE) is a run-level switch measured in
# the smoke test; the choice is frozen in config.yaml per run.

suppressMessages({ library(ANCOMBC); library(TreeSummarizedExperiment); library(S4Vectors) })
source(file.path(Sys.getenv("HARNESS_R_DIR", "."), "lib_common.R"))

WORKDIR <- env_require("WORKDIR")
CELL_ID <- env_require("CELL_ID"); SCENARIO <- env_require("SCENARIO")
REP <- as.integer(env_num("REP")); B <- as.integer(env_num("B"))
BASE_SEED <- env_num("BASE_SEED"); CELL_INDEX <- as.integer(env_num("CELL_INDEX"))
PSEUDO_SENS <- toupper(Sys.getenv("PSEUDO_SENS", "TRUE")) == "TRUE"
ds <- dataset_seed(BASE_SEED, CELL_INDEX, REP)

meta <- read.delim(file.path(WORKDIR, "meta.tsv"))
coldata <- DataFrame(group = factor(meta$group, levels = c("G1", "G2")),
                     row.names = meta$sample_id)
timings <- file.path(WORKDIR, "timings.tsv")

fit_once <- function(mat, psens) {
  tse <- TreeSummarizedExperiment(assays = list(counts = mat), colData = coldata)
  tryCatch(
    ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
             fix_formula = "group", p_adj_method = "BH",
             prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
             struc_zero = FALSE, neg_lb = FALSE, alpha = 0.05,
             global = FALSE, pairwise = FALSE, pseudo_sens = psens,
             n_cl = 1, verbose = FALSE)$res,
    error = function(e) e)
}
# ancombc2 hard-errors when any taxon has zero variance (frequent in rarefied
# tables). Retry excluding the taxa named in the error; they are recorded with
# reason "zero_variance" and count as non-detections (intention-to-test).
# psens: PSEUDO_SENS for raw/rarefied; FALSE inside boots replicates (判断2,
# 2026-08-11 — ANCOMBC 2.12 measures pseudo_sens at ~3x fit cost; the boots
# stability score is itself a robustness metric, so the pseudo-count
# sensitivity filter is kept only where the primary claims live).
fit_one <- function(mat, psens = PSEUDO_SENS) {
  keep <- rownames(mat)[rowSums(mat) > 0]
  zv <- character(0)
  for (attempt in 1:4) {
    use <- setdiff(keep, zv)
    if (!length(use)) return(NULL)
    res <- fit_once(mat[use, , drop = FALSE], psens)
    if (!inherits(res, "condition")) break
    msg <- conditionMessage(res)
    bad <- intersect(use, unlist(strsplit(msg, "[^A-Za-z0-9_.:-]+")))
    if (!length(bad)) {
      cat("ancombc2 error (unparseable):", msg, "\n", file = stderr()); return(NULL)
    }
    zv <- union(zv, bad)
    res <- NULL
  }
  if (is.null(res) || inherits(res, "condition")) {
    cat("ancombc2 error: retries exhausted\n", file = stderr()); return(NULL)
  }
  ss_col <- grep("^passed_ss", colnames(res), value = TRUE)
  out <- data.frame(taxon = res$taxon,
             effect = res[["lfc_groupG2"]], se = res[["se_groupG2"]],
             p = res[["p_groupG2"]], q = res[["q_groupG2"]],
             passed_ss = if (length(ss_col)) as.integer(res[[ss_col[1]]]) else NA_integer_,
             detected = as.integer(!is.na(res[["q_groupG2"]]) & res[["q_groupG2"]] < 0.05),
             stringsAsFactors = FALSE)
  attr(out, "zero_variance") <- zv
  out
}
as_rows <- function(fitdf, taxa, condition) {
  out <- empty_results(taxa, CELL_ID, SCENARIO, REP, "ancombc2", condition)
  if (is.null(fitdf)) { out$reason <- "fit_error"; out$detected <- 0L; return(out) }
  i <- match(taxa, fitdf$taxon)
  out$effect <- fitdf$effect[i]; out$se <- fitdf$se[i]
  out$p <- fitdf$p[i]; out$q <- fitdf$q[i]
  out$passed_ss <- fitdf$passed_ss[i]
  out$detected <- ifelse(is.na(i), 0L, fitdf$detected[i])
  out$reason <- ifelse(is.na(i), "allzero",
                       ifelse(is.na(fitdf$q[i]), "na_q", "ok"))
  out$reason[out$taxon %in% attr(fitdf, "zero_variance")] <- "zero_variance"
  out
}

set.seed(ds + SEED_OFF_METHOD)
rows <- list()
raw <- read_counts(file.path(WORKDIR, "raw.tsv")); taxa <- rownames(raw)

t0 <- proc.time()[3]
rows$raw <- as_rows(fit_one(raw), taxa, "raw")
append_timing(timings, CELL_ID, REP, "ancombc2", "raw", NA, proc.time()[3] - t0)

t0 <- proc.time()[3]
rows$rarefied <- as_rows(fit_one(read_counts(file.path(WORKDIR, "rarefied.tsv"))), taxa, "rarefied")
append_timing(timings, CELL_ID, REP, "ancombc2", "rarefied", NA, proc.time()[3] - t0)

reps_list <- vector("list", B)
for (b in seq_len(B)) {
  t0 <- proc.time()[3]
  f <- fit_one(read_counts(file.path(WORKDIR, "boots", sprintf("b%03d.tsv", b))), psens = FALSE)
  if (!is.null(f)) reps_list[[b]] <- f[, c("taxon", "effect", "q", "detected")]
  append_timing(timings, CELL_ID, REP, "ancombc2", "boots", b, proc.time()[3] - t0)
}
bs <- summarize_boots(reps_list, taxa)
bo <- empty_results(taxa, CELL_ID, SCENARIO, REP, "ancombc2", "boots")
bo$boots_B <- B; bo$boots_ok <- bs$boots_ok
bo$boots_detect_count <- bs$boots_detect_count
bo$boots_q_median <- bs$boots_q_median; bo$boots_q_q25 <- bs$boots_q_q25
bo$boots_q_q75 <- bs$boots_q_q75; bo$boots_effect_median <- bs$boots_effect_median
bo$boots_sign_agree <- bs$boots_sign_agree
bo$effect <- bs$boots_effect_median
bo$detected <- as.integer(bs$boots_ok > 0 & bs$boots_detect_count / pmax(bs$boots_ok, 1L) >= 0.5)
bo$reason <- if (bs$boots_ok[1] == 0) "fit_error" else "ok"
rows$boots <- bo

write_results(do.call(rbind, rows), file.path(WORKDIR, "results_ancombc2.tsv"))
cat(sprintf("ancombc2 done: %s rep %d (pseudo_sens=%s)\n", CELL_ID, REP, PSEUDO_SENS))
