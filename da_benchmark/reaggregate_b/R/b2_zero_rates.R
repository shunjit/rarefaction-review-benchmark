#!/usr/bin/env Rscript
# B-2 (CG-07 後段): realized zero-rate table from per-rep geninfo.
# D0 (2,500 / 10,000 / 40,000) is the pre-registered zero-rate lever; this
# tabulates the realized zero fractions (raw and rarefied tables, tested taxa)
# per cell so the response can cite achieved — not nominal — sparsity.
# usage: Rscript b2_zero_rates.R <dsim_run_dir> <out_dir>

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
run_dir <- args[1]; out_dir <- args[2]
stopifnot(dir.exists(run_dir))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gi <- do.call(rbind, lapply(
  Sys.glob(file.path(run_dir, "results", "*", "rep_[0-9]*_geninfo.tsv")), read.delim))
stopifnot(nrow(gi) > 0)

fmt <- function(x) c(mean = mean(x), sd = sd(x), min = min(x), max = max(x))
grp <- split(gi, gi$cell_id)
out <- do.call(rbind, lapply(grp, function(g) {
  zr <- fmt(g$zero_frac_raw); zf <- fmt(g$zero_frac_rarefied)
  data.frame(cell_id = g$cell_id[1], scenario = g$scenario[1],
    d0 = g$d0[1], rho = g$rho[1], n_per_group = g$n_per_group[1],
    theta_mode = g$theta_mode[1], fc = g$fc[1], da_frac = g$da_frac[1],
    n_reps = nrow(g),
    raref_depth_median = median(g$raref_depth),
    achieved_rho_mean = mean(g$achieved_rho),
    tested_n_mean = mean(g$tested_n),
    zero_raw_mean = zr["mean"], zero_raw_sd = zr["sd"],
    zero_raw_min = zr["min"], zero_raw_max = zr["max"],
    zero_rarefied_mean = zf["mean"], zero_rarefied_sd = zf["sd"],
    zero_rarefied_min = zf["min"], zero_rarefied_max = zf["max"],
    stringsAsFactors = FALSE)
}))
rownames(out) <- NULL
out <- out[order(out$scenario, out$d0, out$cell_id), ]
write.table(out, file.path(out_dir, "b2_realized_zero_rates_by_cell.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# focused view: the D0 lever series (everything else at main-cell settings)
lever <- out[out$cell_id %in% c("null_r100_d2500", "null_r100", "null_r100_d40000",
                                "spike_r100_f3_d2500", "spike_r100_f3",
                                "spike_r100_f3_d40000"), ]
write.table(lever, file.path(out_dir, "b2_realized_zero_rates_d0_lever.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("cells: ", nrow(out), "; reps total: ", nrow(gi))
print(lever[, c("cell_id", "d0", "n_reps", "raref_depth_median",
                "zero_raw_mean", "zero_rarefied_mean")])
