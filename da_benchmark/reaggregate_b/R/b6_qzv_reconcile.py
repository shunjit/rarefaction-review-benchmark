#!/usr/bin/env python3
"""B-6 (CG-17, GM-06): reconcile the phase-1 PERMANOVA qzv artifacts with the
ledger n series and the 14_2_7 primary record, and re-derive the Table S1 (D)
R^2 column from the unrounded pseudo-F (identity R^2 = F(k-1)/[F(k-1)+(n-k)],
k = 4, round-half-up to 3 dp).

The qzv files are the actual QIIME 2 beta-group-significance outputs: sample
size and test statistic are read from <uuid>/data/index.html inside each
artifact, so this is a primary-output check (not a back-calculation).

usage: b6_qzv_reconcile.py <sensitivity_dir> <out_dir>
  <sensitivity_dir> = .../Rice_v2/04_sensitivity_analysis (has depth_*/)
"""
import sys, os, re, html, zipfile, hashlib, decimal

sens_dir, out_dir = sys.argv[1], sys.argv[2]
os.makedirs(out_dir, exist_ok=True)

# reference series being reconciled
LEDGER_N   = {10000: 492, 15000: 490, 20000: 488, 22714: 486}   # sample ledger (基盤 A)
RECORD_F   = {10000: "165.488004", 15000: "165.015866",
              20000: "164.840491", 22714: "165.417978"}          # 14_2_7 primary record
PRINTED_R2 = {10000: "0.504", 15000: "0.504", 20000: "0.505", 22714: "0.507"}  # submitted Table S1 (D)
MISPRINT_N = {10000: 487, 15000: 485, 20000: 483, 22714: 481}    # misprinted series (submitted table)
K = 4

def round_half_up(x, nd):
    return decimal.Decimal(x).quantize(decimal.Decimal("1." + "0" * nd),
                                       rounding=decimal.ROUND_HALF_UP)

rows = []
for depth in sorted(LEDGER_N):
    qzv = os.path.join(sens_dir, f"depth_{depth}", "permanova_wunifrac_compartment.qzv")
    if not os.path.exists(qzv):
        sys.exit(f"ERROR: required qzv MISSING: {qzv} - stopping.")
    with zipfile.ZipFile(qzv) as z:
        name = [n for n in z.namelist() if n.endswith("data/index.html")][0]
        uuid = name.split("/")[0]
        text = re.sub(r"\s+", " ", html.unescape(
            re.sub(r"<[^>]+>", " ", z.read(name).decode("utf-8"))))
    n_qzv  = int(re.search(r"sample size (\d+)", text).group(1))
    f_qzv  = re.search(r"test statistic ([0-9.]+)", text).group(1)
    groups = int(re.search(r"number of groups (\d+)", text).group(1))
    sha = hashlib.sha256(open(qzv, "rb").read()).hexdigest()

    F = decimal.Decimal(f_qzv)
    r2_exact = (F * (K - 1)) / (F * (K - 1) + (n_qzv - K))
    r2_new = round_half_up(r2_exact, 3)
    # the same identity on the misprinted n, for the downgraded side-note
    r2_mis = round_half_up((F * (K - 1)) / (F * (K - 1) + (MISPRINT_N[depth] - K)), 3)

    rows.append(dict(depth=depth, qzv_uuid=uuid, qzv_sha256=sha,
        n_qzv=n_qzv, groups_qzv=groups, pseudo_f_qzv=f_qzv,
        n_ledger=LEDGER_N[depth], n_match=n_qzv == LEDGER_N[depth],
        f_record=RECORD_F[depth], f_match=f_qzv == RECORD_F[depth],
        r2_exact=str(r2_exact.quantize(decimal.Decimal("1.000000"))),
        r2_rederived=str(r2_new), r2_printed=PRINTED_R2[depth],
        r2_printed_match=str(r2_new) == PRINTED_R2[depth],
        n_misprint=MISPRINT_N[depth], r2_on_misprint_n=str(r2_mis)))

hdr = list(rows[0].keys())
with open(os.path.join(out_dir, "b6_qzv_reconciliation.tsv"), "w") as fh:
    fh.write("\t".join(hdr) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[c]) for c in hdr) + "\n")

for r in rows:
    print(f"depth {r['depth']:>6}: n {r['n_qzv']} (ledger {r['n_ledger']}, "
          f"match={r['n_match']}), F {r['pseudo_f_qzv']} (record match={r['f_match']}), "
          f"R2 exact {r['r2_exact']} -> {r['r2_rederived']} "
          f"(printed {r['r2_printed']}, match={r['r2_printed_match']}); "
          f"misprint-n R2 {r['r2_on_misprint_n']}")
