#!/usr/bin/env python3
# Per-genus SGB counts are tested for association with GTDB r220 species-cluster counts.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/sampling_effort_all_genera.py
# Output: results/sampling_effort_all_genera.tsv
import os, sys
from collections import defaultdict

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
MAN = B + "/data/herptile_bacillota_A_HQ_manifest_with_source.tsv"
CDB = B + "/results/drep_herptile_95ani_2229/data_tables/Cdb.csv"
TAX = "/srv/projects/db/gtdbtk/220/taxonomy/bac120_taxonomy_r220_reps.tsv"
PRIOR = B + "/results/genus_expansion.tsv"
OUT = B + "/results/sampling_effort_all_genera.tsv"

FAM = "Ruminococcaceae"
SMALL = 20
NPERM = 9999
SEED = 20260731

for p in (MAN, CDB, TAX):
    if not os.path.exists(p):
        print("MISSING:", p); sys.exit(1)

def fld(tax, pre):
    for f in tax.split(";"):
        f = f.strip()
        if f.startswith(pre):
            return f[len(pre):] if len(f) > len(pre) else "UNASSIGNED"
    return ""

with open(MAN) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
gi = hdr.index("bin_id"); ti = hdr.index("taxonomy")
si = hdr.index("source"); hi = hdr.index("has_metadata")
ci = hdr.index("completeness"); xi = hdr.index("contamination")

rows = []
with open(MAN) as fh:
    fh.readline()
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(gi, ti, si, hi, ci, xi) or not p[gi]:
            continue
        if fld(p[ti], "f__") != FAM:
            continue
        try:
            comp = float(p[ci]); cont = float(p[xi])
        except ValueError:
            comp = cont = float("nan")
        rows.append((p[gi], fld(p[ti], "g__") or "UNASSIGNED",
                     p[si].strip().upper(), p[hi].strip(), comp, cont))
print("%s MAGs in manifest: %d" % (FAM, len(rows)))

def strip_ext(x):
    for e in (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz"):
        if x.endswith(e):
            return x[: -len(e)]
    return x

sgb = {}
with open(CDB) as fh:
    ch = fh.readline().rstrip("\n").replace('"', "").split(",")
    cg = ch.index("genome"); cs = ch.index("secondary_cluster")
    for line in fh:
        p = line.rstrip("\n").replace('"', "").split(",")
        if len(p) > max(cg, cs):
            sgb[strip_ext(p[cg])] = p[cs]

def counts(sel):
    d = defaultdict(set)
    for g, gen, src, hm, comp, cont in sel:
        if g in sgb:
            d[gen].add(sgb[g])
    named = {k: len(v) for k, v in d.items() if k != "UNASSIGNED"}
    return d, named, sum(len(v) for v in d.values()), len(d)

print()
print("=== RECONCILIATION: 276 / 42 vs recorded 274 / 41 ===")
variants = [
    ("no filter", rows),
    ("has_metadata true", [r for r in rows if r[3].lower() in ("true", "yes", "1")]),
    ("completeness >= 90", [r for r in rows if r[4] == r[4] and r[4] >= 90]),
    ("completeness >= 90, contam <= 5",
     [r for r in rows if r[4] == r[4] and r[4] >= 90 and r[5] <= 5]),
    ("completeness > 90", [r for r in rows if r[4] == r[4] and r[4] > 90]),
    ("exclude VIVARIUM", [r for r in rows if r[2] != "VIVARIUM"]),
    ("WILD + ZOO only", [r for r in rows if r[2] in ("WILD", "ZOO")]),
]
print("%-34s %6s %7s %7s %s" % ("variant", "MAGs", "SGBs", "genera", "match"))
hits = []
for lab, sel in variants:
    _, _, t, ng = counts(sel)
    ok = "<== 274/41" if (t == 274 and ng == 41) else ""
    if ok:
        hits.append(lab)
    print("%-34s %6d %7d %7d %s" % (lab, len(sel), t, ng, ok))

if os.path.exists(PRIOR):
    print()
    print("=== PRIOR FILE results/genus_expansion.tsv ===")
    with open(PRIOR) as fh:
        ph = fh.readline().rstrip("\n").split("\t")
        print("columns:", ph)
        n = 0
        for line in fh:
            n += 1
        print("rows:", n)
else:
    print()
    print("results/genus_expansion.tsv NOT FOUND, cannot compare directly")

if hits:
    print()
    print("RECONCILED by: %s" % hits[0])
    lab, sel = [v for v in variants if v[0] == hits[0]][0]
else:
    print()
    print("NO VARIANT REPRODUCES 274/41. Using the unfiltered set and")
    print("reporting 276/42. The 2-SGB gap must be tracked down separately.")
    sel = rows

ref = defaultdict(int)
for line in open(TAX):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 2 and fld(p[1], "f__") == FAM:
        ref[fld(p[1], "g__") or "UNASSIGNED"] += 1
allg = sorted(k for k in ref if k != "UNASSIGNED")
print()
print("GTDB %s: %d genomes across %d named genera" % (sum(ref.values()), len(allg), len(allg)))

import random
rnd = random.Random(SEED)

def spearman(x, y):
    n = len(x)
    def rank(v):
        o = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[o[j + 1]] == v[o[i]]:
                j += 1
            a = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[o[k]] = a
            i = j + 1
        return r
    rx, ry = rank(x), rank(y)
    mx = sum(rx) / n; my = sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx = sum((rx[i] - mx) ** 2 for i in range(n)) ** 0.5
    dy = sum((ry[i] - my) ** 2 for i in range(n)) ** 0.5
    return num / (dx * dy) if dx and dy else float("nan")

def perm_p(x, y, rho):
    ge = 1
    yy = list(y)
    for _ in range(NPERM):
        rnd.shuffle(yy)
        if abs(spearman(x, yy)) >= abs(rho):
            ge += 1
    return float(ge) / (NPERM + 1)

def run(label, subset):
    d, named, tot, ng = counts(subset)
    xa = [ref[g] for g in allg]
    ya = [named.get(g, 0) for g in allg]
    ra = spearman(xa, ya); pa = perm_p(xa, ya, ra)
    obs = sorted(named)
    xo = [ref.get(g, 0) for g in obs]
    yo = [named[g] for g in obs]
    ro = spearman(xo, yo); po = perm_p(xo, yo, ro)
    small = sum(named[g] for g in obs if ref.get(g, 0) <= SMALL)
    den = sum(named.values())
    zero = [g for g in allg if named.get(g, 0) == 0]
    big0 = sorted(zero, key=lambda g: -ref[g])[:8]
    print()
    print("=== %s ===" % label)
    print("  MAGs %d | SGBs %d | genera with >=1 SGB %d" % (len(subset), tot, ng))
    print("  SGBs in genera with <=%d GTDB genomes: %d of %d (%.1f%%)"
          % (SMALL, small, den, 100.0 * small / den if den else 0))
    print("  CORRECTED, all %d GTDB genera incl. zeros : rho = %+.3f, p = %.4f"
          % (len(allg), ra, pa))
    print("  conditioned on recovery (old test), n=%d : rho = %+.3f, p = %.4f"
          % (len(obs), ro, po))
    print("  GTDB genera with ZERO SGBs here: %d" % len(zero))
    print("  largest of them: %s"
          % ", ".join("%s(%d)" % (g, ref[g]) for g in big0))
    return (label, len(subset), tot, ng, small, den, ra, pa, ro, po, len(obs), len(zero))

res = []
res.append(run("ALL ANIMALS", sel))
res.append(run("WILD ONLY", [r for r in sel if r[2] == "WILD"]))

with open(OUT, "w") as f:
    f.write("set\tn_mags\tn_sgbs\tn_genera_with_sgb\tsgbs_in_small_genera\t"
            "sgbs_named\trho_all_genera\tp_all_genera\trho_conditioned\t"
            "p_conditioned\tn_conditioned\tn_zero_genera\tnperm\n")
    for r in res:
        f.write("\t".join(str(x) if not isinstance(x, float) else "%.4f" % x
                          for x in r) + "\t%d\n" % NPERM)
print()
print("wrote", OUT)
print()
print("The corrected test includes every GTDB %s genus, scoring zero where no" % FAM)
print("herptile SGB was recovered. The old test dropped those genera, which")
print("removed exactly the high-GTDB zero-recovery points and inflated rho.")
print("DONE_SAMPLING_EFFORT_ALL")
# SENTINEL_END
