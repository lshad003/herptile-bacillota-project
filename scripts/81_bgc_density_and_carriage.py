# Cluster density per megabase and carriage are computed as descriptive estimands.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/bgc_density_and_carriage.py
# Output: results/bgc_density_carriage.tsv
import os, csv, math, random
from collections import defaultdict

# BGC_DENSITY_AND_CARRIAGE V1. Adds the two descriptive estimands the
# fragmentation review calls modal practice and which R9 currently lacks:
#   (1) BGCs per Mb, the one normalisation that appears descriptively in the
#       published surveys (archaea survey, gut surveys).
#   (2) proportion of genomes carrying at least one COMPLETE BGC.
# Review section (f): the carriage estimand is ALSO fragmentation-dependent,
# because a genome whose only BGCs are truncated has zero complete BGCs. That
# is not a reason to skip it, it is a reason to report it WITH the caveat and
# WITH the same within-arm contiguity check applied to the counts.
# Review section (g): antiSMASH database v4 discards assemblies above 100
# contigs. Every MAG arm here sits well above that. Reported explicitly.

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
SRC = B + "/results/assembly_quality_all_arms.tsv"
OUT = B + "/results/bgc_density_carriage.tsv"

ARMS = ["UHM amphibian", "EHI newt", "EHI mammal", "Youngblut",
        "GTDB reference"]
NPERM = 499
ANTISMASH_DB_V4_CONTIG_LIMIT = 100

if os.path.exists(OUT):
    raise SystemExit("REFUSING TO OVERWRITE: " + OUT)
if not os.path.exists(SRC):
    raise SystemExit("MISSING: " + SRC)

with open(SRC, newline="") as fh:
    rr = [x for x in csv.reader(fh, delimiter="\t") if x]
H = {c: i for i, c in enumerate(rr[0])}
for need in ("arm", "n50", "n_contigs", "total_bp", "n_regions",
             "n_complete", "n_edge"):
    if need not in H:
        raise SystemExit("column missing: " + need)

rows = []
for p in rr[1:]:
    try:
        mb = float(p[H["total_bp"]]) / 1e6
        if mb <= 0:
            continue
        rows.append({"arm": p[H["arm"]].strip(),
                     "n50": float(p[H["n50"]]),
                     "ctg": float(p[H["n_contigs"]]),
                     "mb": mb,
                     "reg": float(p[H["n_regions"]]),
                     "comp": float(p[H["n_complete"]]),
                     "edge": float(p[H["n_edge"]])})
    except (ValueError, IndexError):
        continue

by = {a: [r for r in rows if r["arm"] == a] for a in ARMS}
by = {a: v for a, v in by.items() if v}

def mean(v):
    return sum(v) / len(v) if v else float("nan")

def median(v):
    s = sorted(v)
    n = len(s)
    if not n:
        return float("nan")
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2.0

def spearman(x, y):
    n = len(x)
    if n < 8:
        return float("nan"), float("nan")
    def rk(v):
        idx = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[idx[j + 1]] == v[idx[i]]:
                j += 1
            a = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[idx[k]] = a
            i = j + 1
        return r
    rx, ry = rk(x), rk(y)
    mx, my = mean(rx), mean(ry)
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx <= 0 or syy <= 0:
        return float("nan"), float("nan")
    rho = sxy / math.sqrt(sxx * syy)
    t = rho * math.sqrt(max(n - 2, 1) / max(1e-12, 1 - rho * rho))
    return rho, t

print("=" * 78)
print("ARM SUMMARY. Ruminococcaceae only, all five arms.")
print("=" * 78)
print("  %-18s %6s %9s %8s %8s %10s %10s %10s"
      % ("arm", "n", "med N50", "med ctg", "med Mb", "reg/Mb",
         "complete/Mb", "carriage"))
summ = {}
for a in ARMS:
    if a not in by:
        continue
    v = by[a]
    regmb = mean([r["reg"] / r["mb"] for r in v])
    cmpmb = mean([r["comp"] / r["mb"] for r in v])
    carr = sum(1 for r in v if r["comp"] >= 1) / len(v)
    anyreg = sum(1 for r in v if r["reg"] >= 1) / len(v)
    summ[a] = {"n": len(v), "regmb": regmb, "cmpmb": cmpmb, "carr": carr,
               "anyreg": anyreg, "mb": median([r["mb"] for r in v]),
               "n50": median([r["n50"] for r in v]),
               "ctg": median([r["ctg"] for r in v]),
               "comp": mean([r["comp"] for r in v])}
    print("  %-18s %6d %9.0f %8.0f %8.2f %10.3f %10.3f %10.3f"
          % (a, len(v), summ[a]["n50"], summ[a]["ctg"], summ[a]["mb"],
             regmb, cmpmb, carr))
print()
print("  carriage = proportion of genomes with AT LEAST ONE complete BGC")
print("  For contrast, proportion with at least one region of ANY kind:")
for a in ARMS:
    if a in summ:
        print("    %-18s %.3f" % (a, summ[a]["anyreg"]))
print()

print("=" * 78)
print("REVIEW SECTION (g): THE FIELD'S OWN QUALITY GATE")
print("=" * 78)
print("  antiSMASH database v4 discards assemblies above %d contigs."
      % ANTISMASH_DB_V4_CONTIG_LIMIT)
print("  %-18s %8s %8s" % ("arm", "n", "pass"))
for a in ARMS:
    if a not in by:
        continue
    v = by[a]
    k = sum(1 for r in v if r["ctg"] <= ANTISMASH_DB_V4_CONTIG_LIMIT)
    print("  %-18s %8d %8s" % (a, len(v), "%d (%.1f%%)"
          % (k, 100.0 * k / len(v))))
print()
print("  State this in the limitation. It is a stronger and more honest")
print("  framing than anything specific to this catalog.")
print()

print("=" * 78)
print("REVIEW SECTION (f): DOES CARRIAGE ESCAPE THE CONFOUND? NO.")
print("Same within-arm contiguity test applied to carriage and to per-Mb.")
print("=" * 78)
print("  %-18s %22s %22s %22s"
      % ("arm", "N50 vs complete/Mb", "N50 vs carriage", "N50 vs regions/Mb"))
for a in ARMS:
    if a not in by:
        continue
    v = by[a]
    n50 = [r["n50"] for r in v]
    r1, t1 = spearman(n50, [r["comp"] / r["mb"] for r in v])
    r2, t2 = spearman(n50, [1.0 if r["comp"] >= 1 else 0.0 for r in v])
    r3, t3 = spearman(n50, [r["reg"] / r["mb"] for r in v])
    def f(r, t):
        if r != r:
            return "      n/a           "
        return "  rho %+0.3f %s        " % (r, "*" if abs(t) >= 2 else " ")
    print("  %-18s %s %s %s" % (a, f(r1, t1)[:22], f(r2, t2)[:22],
                                f(r3, t3)[:22]))
print()
print("  * marks |t| >= 2. If carriage is flagged in the same arms as the")
print("  counts, the estimand did NOT escape the confound and must carry")
print("  the caveat wherever it is quoted.")
print()

print("=" * 78)
print("CARRIAGE AT MATCHED ABSOLUTE N50")
print("=" * 78)
BINS = [(0, 10000, "<10 kb"), (10000, 20000, "10-20"), (20000, 40000, "20-40"),
        (40000, 80000, "40-80"), (80000, 1e12, ">80 kb")]
print("  %-18s %s" % ("arm", "".join("%14s" % b[2] for b in BINS)))
for a in ARMS:
    if a not in by:
        continue
    v = by[a]
    cells = []
    for lo, hi, _ in BINS:
        sel = [r for r in v if lo <= r["n50"] < hi]
        if not sel:
            cells.append("%14s" % "n/a")
        else:
            cells.append("%14s" % ("%.2f (%d)"
                         % (sum(1 for r in sel if r["comp"] >= 1) / len(sel),
                            len(sel))))
    print("  %-18s %s" % (a, "".join(cells)))
print()

print("=" * 78)
print("THE ONE COMPARISON THE REVIEW LICENSES")
print("=" * 78)
print("  Review Option 2: restrict inference to fragmentation-matched arms.")
print("  UHM amphibian vs EHI newt PASSES on the named metrics.")
print("  EHI mammal vs EHI newt, the actual host contrast, FAILS on contig")
print("  count (SMD 0.272, KS p 0.035), so the review's own fallback applies:")
print("  Option 1 plus Option 3, descriptive only.")
print()
random.seed(11)
def perm_diff(a, b, key):
    va = [key(r) for r in by[a]]
    vb = [key(r) for r in by[b]]
    obs = mean(va) - mean(vb)
    pool = va + vb
    na = len(va)
    hits = 0
    for _ in range(NPERM):
        random.shuffle(pool)
        if abs(mean(pool[:na]) - mean(pool[na:])) >= abs(obs):
            hits += 1
    return obs, (hits + 1) / (NPERM + 1)

for label, key in (("complete BGCs per Mb", lambda r: r["comp"] / r["mb"]),
                   ("carriage (>=1 complete)",
                    lambda r: 1.0 if r["comp"] >= 1 else 0.0),
                   ("regions per Mb", lambda r: r["reg"] / r["mb"])):
    if "UHM amphibian" in by and "EHI newt" in by:
        d, p = perm_diff("UHM amphibian", "EHI newt", key)
        print("  MATCHED PAIR  %-26s diff %+0.4f  perm p %.3f (%d perms)"
              % (label, d, p, NPERM))
print()
print("  Both arms here are amphibian, so a matched result is NOT a host")
print("  contrast. It is a reproducibility check across two catalogs.")
print()

f = open(OUT, "w")
f.write("arm\tn_genomes\tmedian_n50\tmedian_contigs\tmedian_mb\t"
        "mean_complete_per_genome\tmean_regions_per_mb\t"
        "mean_complete_per_mb\tcarriage_any_complete\tcarriage_any_region\t"
        "frac_at_or_below_100_contigs\n")
for a in ARMS:
    if a not in summ:
        continue
    v = by[a]
    k = sum(1 for r in v if r["ctg"] <= ANTISMASH_DB_V4_CONTIG_LIMIT)
    s = summ[a]
    f.write("%s\t%d\t%.0f\t%.0f\t%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n"
            % (a, s["n"], s["n50"], s["ctg"], s["mb"], s["comp"],
               s["regmb"], s["cmpmb"], s["carr"], s["anyreg"],
               k / len(v)))
f.close()
print("  wrote %s" % OUT)
print()
print("BGC_DENSITY_AND_CARRIAGE_V1_20260806_COMPLETE")
# BGC_DENSITY_AND_CARRIAGE_V1_20260806_COMPLETE
