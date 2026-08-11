# The endotherm arms are added to the contiguity table.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/assembly_quality_endo.py
# Output: results/assembly_quality_all_arms.tsv
# ASSEMBLY_QUALITY_ENDO_V1_20260806
# Adds the 328 endotherm-arm genomes (job 27253466) to the contiguity table
# and asks the question the parse raised: is the Lagomorpha elevation a
# contiguity effect?
#
# WHAT THE PARSE SHOWED: EHI mammal mean complete BGCs 0.56 vs EHI newt 0.33,
# same consortium and pipeline. But by host order:
#   Rodentia   119  0.35     <- on top of the newts
#   Carnivora   82  0.38     <- on top of the newts
#   Lagomorpha  76  1.07     <- three times the others, 27% of the arm
# The arm mean is dominated by one host order. Lagomorpha also has mean ALL
# regions 3.43 vs 2.45 and 1.59, which is consistent with either better
# assemblies or genuinely more BGCs. This script separates those.
#
# METHOD: N50, contig count and length computed from the staged fastas by the
# same code as scripts/assembly_quality_arms.py, so all arms are measured one
# way. A NEW combined file is written; the existing one is not modified.
#
# R9 established that complete-BGC recovery is a function of N50 that does
# not differ by catalog (at matched N50: 0.34/0.14/0.37 in the 20-40 kb bin,
# 0.67/0.65/0.68 in 40-80 kb). If Lagomorpha sits in a higher N50 band, the
# elevation is contiguity. If it is elevated WITHIN an N50 band, it is not.

import os, sys, gzip
import numpy as np

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
ENDO_DIR_LIST = B + "/work/bgc_endo/bgc_endo_input_list.txt"
GEN_ENDO = B + "/results/bgc_per_genome_endo.tsv"
PREV = B + "/results/assembly_quality_arms.tsv"
OUT = B + "/results/assembly_quality_all_arms.tsv"

if os.path.exists(OUT):
    raise SystemExit("REFUSING TO OVERWRITE %s, move it first" % OUT)
for p in (ENDO_DIR_LIST, GEN_ENDO, PREV):
    if not os.path.exists(p):
        sys.exit("MISSING: %s" % p)

try:
    from scipy.stats import ks_2samp
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False


def read_tsv(path):
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                f = f + [""] * (len(header) - len(f))
            rows.append(dict(zip(header, f)))
    return rows


def contig_lengths(path):
    lens = []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                lens.append(0)
            elif lens:
                lens[-1] += len(line.strip())
    return [x for x in lens if x > 0]


def n50(lens):
    if not lens:
        return 0
    s = sorted(lens, reverse=True)
    half = sum(s) / 2.0
    run = 0
    for L in s:
        run += L
        if run >= half:
            return L
    return s[-1]


paths = {}
with open(ENDO_DIR_LIST) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            paths[f[1].strip()] = f[0].strip()
print("endotherm input list: %d entries" % len(paths))

qual = {}
for i, (safe, p) in enumerate(sorted(paths.items())):
    if not os.path.exists(p):
        print("MISSING FASTA: %s" % p)
        continue
    lens = contig_lengths(p)
    qual[safe] = dict(n_contigs=len(lens), total_bp=sum(lens), n50=n50(lens),
                      longest=max(lens) if lens else 0)
    if (i + 1) % 100 == 0:
        print("   %d measured" % (i + 1))
print("endotherm genomes measured: %d" % len(qual))

endo = read_tsv(GEN_ENDO)
ARMLAB = {"ehi_mammal": "EHI mammal", "youngblut": "Youngblut"}
rows = []
for r in endo:
    q = qual.get(r["genome"])
    if q is None:
        continue
    arm = ARMLAB.get(r["arm"], r["arm"])
    if arm == "EHI mammal" and r.get("host_order", ""):
        sub = "EHI mammal: " + r["host_order"]
    else:
        sub = arm
    rows.append(dict(genome=r["genome"], arm=arm, subgroup=sub,
                     genus=r.get("genus", ""),
                     host_order=r.get("host_order", ""),
                     n_regions=int(r["n_regions"]),
                     n_complete=int(r["n_complete"]),
                     n_edge=int(r["n_edge"]), **q))

prev = read_tsv(PREV)
for r in prev:
    rows.append(dict(genome=r["genome"], arm=r["arm"], subgroup=r["arm"],
                     genus=r.get("genus", ""), host_order="",
                     n_contigs=int(r["n_contigs"]), total_bp=int(r["total_bp"]),
                     n50=int(r["n50"]), longest=int(r["longest_contig"]),
                     n_regions=int(r["n_regions"]),
                     n_complete=int(r["n_complete"]),
                     n_edge=int(r["n_edge"])))
print("combined rows: %d" % len(rows))

ORDER = ["UHM amphibian", "EHI newt", "EHI mammal", "Youngblut",
         "GTDB reference"]
byarm = {a: [r for r in rows if r["arm"] == a] for a in ORDER}

print()
print("=" * 78)
print("ASSEMBLY QUALITY, ALL FIVE ARMS, Ruminococcaceae only")
print("=" * 78)
print("  %-16s %6s %11s %12s %10s %14s"
      % ("arm", "n", "med N50", "med contigs", "med Mb", "mean COMPLETE"))
for a in ORDER:
    v = byarm[a]
    if not v:
        continue
    print("  %-16s %6d %11.0f %12.0f %10.2f %14.2f"
          % (a, len(v), np.median([r["n50"] for r in v]),
             np.median([r["n_contigs"] for r in v]),
             np.median([r["total_bp"] for r in v]) / 1e6,
             np.mean([r["n_complete"] for r in v])))

print()
print("=" * 78)
print("THE KEY CONTRAST: EHI mammal vs EHI newt, matched pipeline")
print("=" * 78)


def smd(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    sp = np.sqrt((x.var(ddof=1) + y.var(ddof=1)) / 2.0)
    return abs(x.mean() - y.mean()) / sp if sp > 0 else 0.0


va, vb = byarm["EHI mammal"], byarm["EHI newt"]
for key, lab in (("n50", "N50"), ("n_contigs", "contig count"),
                 ("total_bp", "assembly length")):
    x = [r[key] for r in va]
    y = [r[key] for r in vb]
    s = smd(x, y)
    verdict = ("well balanced" if s < 0.1 else
               "tolerable" if s < 0.25 else "NOT BALANCED")
    line = "  %-18s SMD %6.3f   %s" % (lab, s, verdict)
    if HAVE_SCIPY:
        st, p = ks_2samp(x, y)
        line += "   KS D=%.3f p=%.2g" % (st, p)
    print(line)

print()
print("=" * 78)
print("IS THE LAGOMORPHA ELEVATION A CONTIGUITY EFFECT?")
print("=" * 78)
subs = ["EHI mammal: Rodentia", "EHI mammal: Carnivora",
        "EHI mammal: Lagomorpha", "EHI newt"]
print("  %-24s %6s %11s %12s %14s %12s"
      % ("group", "n", "med N50", "med contigs", "mean COMPLETE", "mean all"))
for s in subs:
    v = [r for r in rows if r["subgroup"] == s]
    if len(v) < 5:
        continue
    print("  %-24s %6d %11.0f %12.0f %14.2f %12.2f"
          % (s, len(v), np.median([r["n50"] for r in v]),
             np.median([r["n_contigs"] for r in v]),
             np.mean([r["n_complete"] for r in v]),
             np.mean([r["n_regions"] for r in v])))

print()
print("  SAME ABSOLUTE N50 BINS, mean complete BGCs (n in brackets).")
print("  This is the test: if Lagomorpha is elevated WITHIN a bin, contiguity")
print("  does not explain it. If it only sits in higher bins, it does.")
edges = [0, 10000, 20000, 40000, 80000, 1e12]
lab = ["<10k", "10-20k", "20-40k", "40-80k", ">80k"]
print("  %-24s %s" % ("group", " ".join("%11s" % x for x in lab)))
for s in subs + ["UHM amphibian", "GTDB reference"]:
    v = [r for r in rows if r["subgroup"] == s or r["arm"] == s]
    if len(v) < 5:
        continue
    out = []
    for j in range(5):
        b = [r["n_complete"] for r in v
             if edges[j] <= r["n50"] < edges[j + 1]]
        out.append("%5.2f(%3d)" % (np.mean(b), len(b)) if len(b) >= 5
                   else "     -     ")
    print("  %-24s %s" % (s, " ".join("%11s" % x for x in out)))

print()
print("  CONTIG-EDGE RATE BY GROUP:")
for s in subs:
    v = [r for r in rows if r["subgroup"] == s]
    if len(v) < 5:
        continue
    tot = sum(r["n_regions"] for r in v)
    e = sum(r["n_edge"] for r in v)
    print("    %-24s %5.1f%%  (%d of %d regions)"
          % (s, 100.0 * e / tot if tot else 0.0, e, tot))

with open(OUT, "w") as f:
    f.write("genome\tarm\tsubgroup\tgenus\thost_order\tn_contigs\ttotal_bp\t"
            "n50\tlongest_contig\tn_regions\tn_complete\tn_edge\n")
    for r in sorted(rows, key=lambda z: (z["arm"], z["genome"])):
        f.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
                % (r["genome"], r["arm"], r["subgroup"], r["genus"],
                   r["host_order"], r["n_contigs"], r["total_bp"], r["n50"],
                   r["longest"], r["n_regions"], r["n_complete"], r["n_edge"]))
print()
print("wrote %s" % OUT)
print("ASSEMBLY_QUALITY_ENDO_V1_20260806_COMPLETE")
