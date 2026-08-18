#!/usr/bin/env python3
"""A-3 (CG-19): extract PERMANOVA / alpha Kruskal-Wallis results from the
common-486 re-run qzv artifacts AND from the phase-1 sensitivity qzv artifacts,
into one long-format table per family, with the R^2 identity re-derivation
(R^2 = F(k-1)/[F(k-1)+(n-k)], k = 4, round-half-up 3 dp) used by B-6.

The qzv files are QIIME 2 beta-group-significance / alpha-group-significance
outputs; sample size, test statistic and p-value are read from
<uuid>/data/index.html inside each artifact (primary output, no back-calc).

usage: extract_a3_permanova.py <run_output_dir> <summary_dir>
  <run_output_dir> = .../a3_common486_v1/output  (has depth_*/)
  phase-1 dir is resolved from $RAW_ROOT (default /Volumes/PS3000/benchmark_data)
"""
import sys, os, re, html, zipfile, decimal

run_out, summary_dir = sys.argv[1], sys.argv[2]
RAW_ROOT = os.environ.get("RAW_ROOT", "/Volumes/PS3000/benchmark_data")
PH1 = os.path.join(RAW_ROOT, "Rice_v2", "04_sensitivity_analysis")
os.makedirs(summary_dir, exist_ok=True)

DEPTHS = [10000, 15000, 20000, 22714]
K = 4

BETA_FILES = {"weighted_unifrac": "permanova_wunifrac_compartment.qzv",
              "unweighted_unifrac": "permanova_uunifrac_compartment.qzv"}
ALPHA_FILES = {"shannon": "shannon_group_significance.qzv",
               "faith_pd": "faith_pd_group_significance.qzv"}


def qzv_text(path):
    with zipfile.ZipFile(path) as z:
        name = [n for n in z.namelist() if n.endswith("data/index.html")][0]
        return re.sub(r"\s+", " ", html.unescape(
            re.sub(r"<[^>]+>", " ", z.read(name).decode("utf-8"))))


def round_half_up(x, nd):
    return decimal.Decimal(x).quantize(decimal.Decimal("1." + "0" * nd),
                                       rounding=decimal.ROUND_HALF_UP)


def parse_beta(path):
    t = qzv_text(path)
    n = int(re.search(r"sample size (\d+)", t).group(1))
    g = int(re.search(r"number of groups (\d+)", t).group(1))
    f = re.search(r"test statistic ([0-9.]+)", t).group(1)
    p = re.search(r"p-value ([0-9.eE-]+)", t).group(1)
    F = decimal.Decimal(f)
    r2 = round_half_up((F * (K - 1)) / (F * (K - 1) + (n - K)), 3)
    return dict(n=n, groups=g, pseudo_F=f, p_value=p, r2_identity=str(r2))


def parse_alpha(path):
    # alpha-group-significance renders the all-groups Kruskal-Wallis result from
    # data/column-<col>.jsonp (index.html cells are empty placeholders), so read
    # the kwAll payload {"H": ..., "p": ...} and the {"initial","filtered"} n.
    with zipfile.ZipFile(path) as z:
        name = [n for n in z.namelist()
                if n.endswith("data/column-compartment.jsonp")][0]
        txt = z.read(name).decode("utf-8")
    h = re.search(r'"H":\s*([0-9.eE+-]+)', txt)
    p = re.search(r'"H":[^}]*"p":\s*([0-9.eE+-]+)', txt)
    n = re.search(r'"filtered":\s*(\d+)', txt)
    return dict(n=n.group(1) if n else "NA",
                H=h.group(1) if h else "NA",
                p_value=p.group(1) if p else "NA")


beta_rows, alpha_rows = [], []
for depth in DEPTHS:
    for series, base in (("common486", os.path.join(run_out, f"depth_{depth}")),
                         ("phase1", os.path.join(PH1, f"depth_{depth}"))):
        for metric, fname in BETA_FILES.items():
            path = os.path.join(base, fname)
            if not os.path.exists(path):
                if series == "common486":
                    sys.exit(f"ERROR: required qzv MISSING: {path}")
                continue  # phase-1 uu may be absent for some depth; record what exists
            beta_rows.append(dict(metric=metric, depth=depth, series=series,
                                  **parse_beta(path)))
        for metric, fname in ALPHA_FILES.items():
            path = os.path.join(base, fname)
            if not os.path.exists(path):
                continue
            alpha_rows.append(dict(metric=metric, depth=depth, series=series,
                                   **parse_alpha(path)))

for name, rows in (("a3_permanova_summary.tsv", beta_rows),
                   ("a3_alpha_summary.tsv", alpha_rows)):
    hdr = list(rows[0].keys())
    with open(os.path.join(summary_dir, name), "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in hdr) + "\n")

for r in beta_rows:
    print(f"{r['metric']:>20} d={r['depth']:>6} {r['series']:>10}: "
          f"n={r['n']} F={r['pseudo_F']} p={r['p_value']} R2={r['r2_identity']}")
