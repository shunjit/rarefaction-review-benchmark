# Reality check — Rice-calibrated DM generator vs Rhizoplane top-300 (both at depth 5,000)

- parameters: pi = empirical mean composition (top-300, coverage 0.441 of reads); theta = 52.7, calibrated to match the real zero fraction (DM MLE reference: 86.3 — shrinks rare taxa and understates sparsity, kept for the record)
- samples compared: 132 (real samples with >= 5000 reads on top-300 ASVs)
- zero fraction: real 0.564 vs sim 0.563  [tuned statistic]
- Shannon (Q1/median/Q3): real 3.525 / 3.985 / 4.297 vs sim 3.732 / 3.853 / 3.923  [untuned]
- prevalence, rank-matched real-vs-sim Spearman: 0.449  [untuned]
- log10 mean count, real-vs-sim Pearson: 0.921  [near-construction: pi is the empirical mean]
- mean-variance: real median var/mean 52.8 vs sim 81.5 (taxa with mean >= 1)  [untuned]

Known limitations (recorded for the response letter): a single-component DM cannot
reproduce the full between-sample heterogeneity of Rhizoplane (field vs greenhouse
runs), so the simulated Shannon spread is narrower than the real one, taxon-level
prevalence ranks decouple from mean-abundance ranks in real data but not under DM
(hence the modest prevalence Spearman), and matching the zero fraction overshoots
the real var/mean ratio by ~1.5x. Mean abundance and global sparsity were
prioritised because effect sizes and rarefaction dropout act on those axes; the
model family (Dirichlet + multinomial) is the one specified by Reviewer #2.

Per-taxon and per-sample tables: reality_per_taxon.tsv / reality_per_sample.tsv.
Generator: p ~ Dirichlet(theta * pi), x ~ Multinomial(5000, p); seeds 20260811/20260812.
