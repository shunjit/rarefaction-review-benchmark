#!/usr/bin/env Rscript
# B-5 (CG-16): attainment of C >= 0.99 at the ADOPTED depths, overall and by
# group, from the frozen per-sample coverage TSVs (基盤 C 出力). No new
# computation — E[C(d)] columns at the adopted depths already exist.
#   Rice 20,000 (ecov_d20000, group = compartment)
#   Soil 20,000 (ecov_d20000, group = treatment)
#   Bioethanol 35,000 (ecov_d35000, group = batch)
# Samples with table_depth < adopted depth have NA ecov (dropped at that
# depth) and are excluded from the attainment denominator but counted.
# usage: Rscript b5_coverage_attainment.R <coverage_output_dir> <out_dir>

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
cov_dir <- args[1]; out_dir <- args[2]
stopifnot(dir.exists(cov_dir))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

spec <- list(
  Rice       = list(file = "rice_coverage_per_sample.tsv",
                    depth = 20000, ecol = "ecov_d20000", group = "compartment"),
  Soil       = list(file = "soil_coverage_per_sample.tsv",
                    depth = 20000, ecol = "ecov_d20000", group = "treatment"),
  Bioethanol = list(file = "bioethanol_coverage_per_sample.tsv",
                    depth = 35000, ecol = "ecov_d35000", group = "batch"))

q1 <- function(x) quantile(x, 0.25, names = FALSE)
q3 <- function(x) quantile(x, 0.75, names = FALSE)
summ_row <- function(ds, depth, grp_label, e) {
  data.frame(dataset = ds, adopted_depth = depth, group = grp_label,
    n = length(e), n_ge_099 = sum(e >= 0.99),
    frac_ge_099 = mean(e >= 0.99),
    ecov_min = min(e), ecov_q25 = q1(e), ecov_median = median(e),
    ecov_q75 = q3(e), ecov_max = max(e), stringsAsFactors = FALSE)
}

overall <- list(); by_grp <- list()
for (ds in names(spec)) {
  s <- spec[[ds]]
  d <- read.delim(file.path(cov_dir, s$file), check.names = FALSE,
                  stringsAsFactors = FALSE)
  stopifnot(s$ecol %in% names(d), s$group %in% names(d))
  e_all <- d[[s$ecol]]
  retained <- !is.na(e_all)
  stopifnot(all(d$table_depth[!retained] < s$depth))   # NA iff dropped
  o <- summ_row(ds, s$depth, "(all)", e_all[retained])
  o$n_total <- nrow(d); o$n_dropped_at_depth <- sum(!retained)
  overall[[ds]] <- o
  g <- d[[s$group]][retained]
  by_grp[[ds]] <- do.call(rbind, lapply(split(e_all[retained], g), function(e2)
    summ_row(ds, s$depth, NA, e2)))
  by_grp[[ds]]$group <- names(split(e_all[retained], g))
}
ov <- do.call(rbind, overall); rownames(ov) <- NULL
bg <- do.call(rbind, by_grp);  rownames(bg) <- NULL
write.table(ov, file.path(out_dir, "b5_c99_attainment_by_dataset.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bg, file.path(out_dir, "b5_c99_attainment_by_group.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(ov[, c("dataset", "adopted_depth", "n", "n_dropped_at_depth",
             "frac_ge_099", "ecov_min", "ecov_median")])
print(bg[, c("dataset", "group", "n", "frac_ge_099", "ecov_min", "ecov_median")])
