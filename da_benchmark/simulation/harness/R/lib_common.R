# Shared helpers for the D-sim harness (D_HARNESS_DESIGN.md v2).
# Sourced by generate_rep.R / fit_*.R. Base R only.

env_require <- function(name) {
  v <- Sys.getenv(name, unset = "")
  if (!nzchar(v)) { cat(sprintf("ERROR: env %s MISSING - stopping.\n", name), file = stderr()); quit(status = 1) }
  v
}
env_num <- function(name) {
  v <- suppressWarnings(as.numeric(env_require(name)))
  if (is.na(v)) { cat(sprintf("ERROR: env %s is not numeric - stopping.\n", name), file = stderr()); quit(status = 1) }
  v
}

# Hierarchical seeds (design section 12)
dataset_seed <- function(base_seed, cell_index, rep) {
  s <- base_seed + cell_index * 1e7 + rep * 1e3
  stopifnot(s < .Machine$integer.max)
  as.integer(s)
}
SEED_OFF_SINGLE_RAREF <- 901L
SEED_OFF_PERMUTE     <- 902L
SEED_OFF_METHOD      <- 903L

read_counts <- function(path) {
  m <- as.matrix(read.delim(path, row.names = 1, check.names = FALSE))
  storage.mode(m) <- "integer"
  m
}
write_counts <- function(m, path) {
  df <- data.frame(taxon = rownames(m), m, check.names = FALSE)
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Exact rarefaction (multivariate hypergeometric, sequential): same law as
# sampling `depth` reads without replacement from the observed reads.
rarefy_vec <- function(x, depth) {
  stopifnot(sum(x) >= depth)
  out <- integer(length(x)); n_rem <- depth; N_rem <- sum(x)
  for (j in seq_along(x)) {
    if (n_rem == 0L) break
    h <- rhyper(1, x[j], N_rem - x[j], n_rem)
    out[j] <- h; n_rem <- n_rem - h; N_rem <- N_rem - x[j]
  }
  out
}
rarefy_mat <- function(m, depth) {
  r <- vapply(seq_len(ncol(m)), function(j) rarefy_vec(m[, j], depth), integer(nrow(m)))
  dimnames(r) <- dimnames(m)
  r
}

# Unified per-taxon result schema (one row per taxon x method x condition)
RESULT_COLS <- c("cell_id", "scenario", "rep", "method", "condition", "taxon",
                 "effect", "se", "p", "q", "detected", "reason", "passed_ss",
                 "boots_B", "boots_ok", "boots_detect_count",
                 "boots_q_median", "boots_q_q25", "boots_q_q75",
                 "boots_effect_median", "boots_sign_agree")

empty_results <- function(taxa, cell_id, scenario, rep, method, condition) {
  df <- data.frame(cell_id = cell_id, scenario = scenario, rep = rep,
                   method = method, condition = condition, taxon = taxa,
                   effect = NA_real_, se = NA_real_, p = NA_real_, q = NA_real_,
                   detected = NA_integer_, reason = "ok", passed_ss = NA_integer_,
                   boots_B = NA_integer_, boots_ok = NA_integer_,
                   boots_detect_count = NA_integer_,
                   boots_q_median = NA_real_, boots_q_q25 = NA_real_, boots_q_q75 = NA_real_,
                   boots_effect_median = NA_real_, boots_sign_agree = NA_real_,
                   stringsAsFactors = FALSE)
  df[, RESULT_COLS]
}

write_results <- function(df, path) {
  stopifnot(identical(colnames(df), RESULT_COLS))
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

append_timing <- function(path, cell_id, rep, method, condition, b_index, elapsed) {
  row <- sprintf("%s\t%s\t%s\t%s\t%s\t%.3f", cell_id, rep, method, condition,
                 ifelse(is.na(b_index), "NA", b_index), elapsed)
  if (!file.exists(path))
    writeLines("cell_id\trep\tmethod\tcondition\tb_index\telapsed_sec", path)
  cat(row, "\n", sep = "", file = path, append = TRUE)
}

# Aggregate B per-replicate fits into one 'boots' row per taxon.
# reps_list: list over replicates of data.frame(taxon, effect, q, detected(0/1))
summarize_boots <- function(reps_list, taxa) {
  ok <- Filter(Negate(is.null), reps_list)
  B_ok <- length(ok)
  out <- data.frame(taxon = taxa, boots_detect_count = 0L, boots_ok = B_ok,
                    boots_q_median = NA_real_, boots_q_q25 = NA_real_,
                    boots_q_q75 = NA_real_, boots_effect_median = NA_real_,
                    boots_sign_agree = NA_real_, stringsAsFactors = FALSE)
  if (B_ok == 0L) return(out)
  qm <- sapply(ok, function(d) d$q[match(taxa, d$taxon)])
  em <- sapply(ok, function(d) d$effect[match(taxa, d$taxon)])
  dm <- sapply(ok, function(d) d$detected[match(taxa, d$taxon)])
  out$boots_detect_count <- as.integer(rowSums(dm == 1L, na.rm = TRUE))
  out$boots_q_median <- apply(qm, 1, median, na.rm = TRUE)
  out$boots_q_q25 <- apply(qm, 1, quantile, probs = 0.25, na.rm = TRUE, names = FALSE)
  out$boots_q_q75 <- apply(qm, 1, quantile, probs = 0.75, na.rm = TRUE, names = FALSE)
  out$boots_effect_median <- apply(em, 1, median, na.rm = TRUE)
  med_sign <- sign(out$boots_effect_median)
  out$boots_sign_agree <- vapply(seq_along(taxa), function(i) {
    e <- em[i, ]; e <- e[!is.na(e)]
    if (!length(e) || med_sign[i] == 0) return(NA_real_)
    mean(sign(e) == med_sign[i])
  }, numeric(1))
  out
}
