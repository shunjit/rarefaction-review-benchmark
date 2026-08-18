#!/usr/bin/env Rscript
# MaAsLin3 runner (D_HARNESS_DESIGN.md v2 section 11).
# Fixed config: formula "~ group", normalization TSS, transform LOG,
# correction BH, min_abundance/min_prevalence/min_variance 0,
# median_comparison_abundance TRUE (documented default), plots off, cores 1.
# Primary call: abundance-model qval_individual < 0.05 (method "maaslin3").
# Prevalence-model rows are stored as method "maaslin3_prev" (secondary).

suppressMessages(library(maaslin3))
source(file.path(Sys.getenv("HARNESS_R_DIR", "."), "lib_common.R"))

WORKDIR <- env_require("WORKDIR")
CELL_ID <- env_require("CELL_ID"); SCENARIO <- env_require("SCENARIO")
REP <- as.integer(env_num("REP")); B <- as.integer(env_num("B"))
BASE_SEED <- env_num("BASE_SEED"); CELL_INDEX <- as.integer(env_num("CELL_INDEX"))
ds <- dataset_seed(BASE_SEED, CELL_INDEX, REP)

meta <- read.delim(file.path(WORKDIR, "meta.tsv"))
metadata <- data.frame(group = factor(meta$group, levels = c("G1", "G2")),
                       row.names = meta$sample_id)
timings <- file.path(WORKDIR, "timings.tsv")

fit_one <- function(mat) {
  keep <- rowSums(mat) > 0
  if (!any(keep)) return(NULL)
  input <- as.data.frame(t(mat[keep, , drop = FALSE]))
  outdir <- file.path(WORKDIR, "m3_tmp")
  unlink(outdir, recursive = TRUE)
  ok <- tryCatch({
    invisible(capture.output(suppressMessages(suppressWarnings(
      maaslin3(input_data = input, input_metadata = metadata, output = outdir,
               formula = "~ group", normalization = "TSS", transform = "LOG",
               correction = "BH", min_abundance = 0, min_prevalence = 0,
               min_variance = 0, median_comparison_abundance = TRUE,
               standardize = TRUE, plot_summary_plot = FALSE,
               plot_associations = FALSE, save_models = FALSE,
               cores = 1, verbosity = "WARN")))))
    TRUE
  }, error = function(e) { cat("maaslin3 error:", conditionMessage(e), "\n", file = stderr()); FALSE })
  if (!ok) return(NULL)
  res_path <- file.path(outdir, "all_results.tsv")
  if (!file.exists(res_path)) return(NULL)
  res <- read.delim(res_path)
  res <- res[res$metadata == "group", , drop = FALSE]
  err <- if ("error" %in% colnames(res)) !is.na(res$error) & nzchar(res$error) else rep(FALSE, nrow(res))
  data.frame(taxon = res$feature, model = res$model,
             effect = res$coef, se = res$stderr,
             p = res$pval_individual, q = res$qval_individual,
             fit_err = err,
             detected = as.integer(!err & !is.na(res$qval_individual) & res$qval_individual < 0.05),
             stringsAsFactors = FALSE)
}
as_rows <- function(fitdf, taxa, condition, model_name, method_label) {
  out <- empty_results(taxa, CELL_ID, SCENARIO, REP, method_label, condition)
  if (is.null(fitdf)) { out$reason <- "fit_error"; out$detected <- 0L; return(out) }
  sub <- fitdf[fitdf$model == model_name, , drop = FALSE]
  i <- match(taxa, sub$taxon)
  out$effect <- sub$effect[i]; out$se <- sub$se[i]
  out$p <- sub$p[i]; out$q <- sub$q[i]
  out$detected <- ifelse(is.na(i), 0L, sub$detected[i])
  out$reason <- ifelse(is.na(i), "allzero",
                       ifelse(sub$fit_err[i], "fit_error",
                              ifelse(is.na(sub$q[i]), "na_q", "ok")))
  out
}

set.seed(ds + SEED_OFF_METHOD)
rows <- list()
raw <- read_counts(file.path(WORKDIR, "raw.tsv")); taxa <- rownames(raw)

t0 <- proc.time()[3]
f <- fit_one(raw)
rows$raw <- as_rows(f, taxa, "raw", "abundance", "maaslin3")
rows$raw_prev <- as_rows(f, taxa, "raw", "prevalence", "maaslin3_prev")
append_timing(timings, CELL_ID, REP, "maaslin3", "raw", NA, proc.time()[3] - t0)

t0 <- proc.time()[3]
f <- fit_one(read_counts(file.path(WORKDIR, "rarefied.tsv")))
rows$rarefied <- as_rows(f, taxa, "rarefied", "abundance", "maaslin3")
rows$rarefied_prev <- as_rows(f, taxa, "rarefied", "prevalence", "maaslin3_prev")
append_timing(timings, CELL_ID, REP, "maaslin3", "rarefied", NA, proc.time()[3] - t0)

reps_list <- vector("list", B)
for (b in seq_len(B)) {
  t0 <- proc.time()[3]
  f <- fit_one(read_counts(file.path(WORKDIR, "boots", sprintf("b%03d.tsv", b))))
  if (!is.null(f)) {
    ab <- f[f$model == "abundance", c("taxon", "effect", "q", "detected")]
    reps_list[[b]] <- ab
  }
  append_timing(timings, CELL_ID, REP, "maaslin3", "boots", b, proc.time()[3] - t0)
}
bs <- summarize_boots(reps_list, taxa)
bo <- empty_results(taxa, CELL_ID, SCENARIO, REP, "maaslin3", "boots")
bo$boots_B <- B; bo$boots_ok <- bs$boots_ok
bo$boots_detect_count <- bs$boots_detect_count
bo$boots_q_median <- bs$boots_q_median; bo$boots_q_q25 <- bs$boots_q_q25
bo$boots_q_q75 <- bs$boots_q_q75; bo$boots_effect_median <- bs$boots_effect_median
bo$boots_sign_agree <- bs$boots_sign_agree
bo$effect <- bs$boots_effect_median
bo$detected <- as.integer(bs$boots_ok > 0 & bs$boots_detect_count / pmax(bs$boots_ok, 1L) >= 0.5)
bo$reason <- if (bs$boots_ok[1] == 0) "fit_error" else "ok"
rows$boots <- bo

unlink(file.path(WORKDIR, "m3_tmp"), recursive = TRUE)
write_results(do.call(rbind, rows), file.path(WORKDIR, "results_maaslin3.tsv"))
cat(sprintf("maaslin3 done: %s rep %d\n", CELL_ID, REP))
