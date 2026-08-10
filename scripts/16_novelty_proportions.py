#!/usr/bin/env python3
# Novelty proportions against GTDB r220 and per-genus expansion are computed for the catalog.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/novelty_proportions.py
# Output: results/novelty_proportions.tsv, results/genus_expansion.tsv
import os, sys
from collections import Counter
import numpy as np

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
SGB = B + "/data/sgb_manifest.tsv"
TAX = "/srv/projects/db/gtdbtk/220/taxonomy/bac120_taxonomy_r220_reps.tsv"
OUT = B + "/results/novelty_proportions.tsv"
OUTG = B + "/results/genus_expansion.tsv"

FAM = "Ruminococcaceae"
UNK = ("UNASSIGNED", "(unassigned)", "")
NPERM = 9999
SEED = 20260803

if not os.path.exists(SGB):
    print("MISSING:", SGB); sys.exit(1)

rows = []
with open(SGB) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    I = {k: i for i, k in enumerate(h)}
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(I.values()):
            continue
        rows.append(dict(sgb=p[I["sgb"]], fam=p[I["family"]].strip(),
                         gen=p[I["genus"]].strip(), sp=p[I["species"]].strip(),
                         wild=p[I["has_wild"]] == "yes",
                         n=int(p[I["n_mags"]])))
print("SGBs in manifest: %d" % len(rows))

def novel_sp(r):
    return r["sp"] in UNK
def novel_gen(r):
    return r["gen"] in UNK
def novel_fam(r):
    return r["fam"] in UNK

def block(sel, label):
    if not sel:
        print("\n%s: empty" % label); return None
    n = len(sel)
    ns = sum(1 for r in sel if novel_sp(r))
    ng = sum(1 for r in sel if novel_gen(r))
    nf = sum(1 for r in sel if novel_fam(r))
    print()
    print("=== %s ===" % label)
    print("  SGBs                              : %d" % n)
    print("  novel at species level (no s__)   : %d (%.1f%%)" % (ns, 100.0*ns/n))
    print("  novel at genus level   (no g__)   : %d (%.1f%%)" % (ng, 100.0*ng/n))
    print("  novel at family level  (no f__)   : %d (%.1f%%)" % (nf, 100.0*nf/n))
    print("  matched an existing GTDB species  : %d (%.1f%%)"
          % (n - ns, 100.0*(n-ns)/n))
    return (label, n, ns, ng, nf)

res = []
res.append(block(rows, "ALL SGBs, all animals"))
res.append(block([r for r in rows if r["wild"]], "ALL SGBs, wild only"))
rs = [r for r in rows if r["fam"] == FAM]
res.append(block(rs, "%s SGBs, all animals" % FAM))
rw = [r for r in rs if r["wild"]]
res.append(block(rw, "%s SGBs, wild only" % FAM))

print()
print("=== BY FAMILY, wild only, top 10 ===")
wf = [r for r in rows if r["wild"]]
print("  %-26s %5s %8s %8s" % ("family", "SGBs", "novel_sp", "%"))
for f, c in Counter(r["fam"] for r in wf).most_common(10):
    sel = [r for r in wf if r["fam"] == f]
    ns = sum(1 for r in sel if novel_sp(r))
    print("  %-26s %5d %8d %7.1f%%"
          % (f if f not in UNK else "(no family)", c, ns, 100.0*ns/c))

with open(OUT, "w") as f:
    f.write("set\tn_sgbs\tnovel_species\tnovel_genus\tnovel_family\t"
            "pct_novel_species\tpct_novel_genus\tpct_novel_family\n")
    for r in res:
        if r is None:
            continue
        lab, n, ns, ng, nf = r
        f.write("%s\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\n"
                % (lab, n, ns, ng, nf, 100.0*ns/n, 100.0*ng/n, 100.0*nf/n))
print()
print("wrote", OUT)

def fld(t, pre):
    for z in t.split(";"):
        z = z.strip()
        if z.startswith(pre):
            return z[len(pre):]
    return ""

ref = Counter()
with open(TAX) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2 and fld(p[1], "f__") == FAM:
            g = fld(p[1], "g__")
            if g:
                ref[g] += 1

wg = Counter(r["gen"] for r in rw if r["gen"] not in UNK)
tot = sum(wg.values())
print()
print("=== GENUS EXPANSION, wild %s SGBs vs GTDB r220 species clusters ===" % FAM)
print("  %-22s %6s %8s %8s" % ("genus", "SGBs", "GTDB", "ratio"))
grows = []
for g, c in wg.most_common():
    nr = ref.get(g, 0)
    ratio = float(c) / nr if nr else float("inf")
    grows.append((g, c, nr, ratio))
    print("  %-22s %6d %8d %8s"
          % (g, c, nr, "%.2f" % ratio if nr else "inf"))

rng = np.random.RandomState(SEED)
gen_list = sorted(ref)
w = np.array([ref[g] for g in gen_list], dtype=float)
w = w / w.sum()
obs_hi = sum(1 for g, c, nr, x in grows if nr > 0 and x >= 2.0)
obs_max = max((x for g, c, nr, x in grows if nr > 0), default=0.0)
cnt_hi = 1; cnt_max = 1
for _ in range(NPERM):
    draw = rng.multinomial(tot, w)
    rr = [draw[i] / float(ref[gen_list[i]]) for i in range(len(gen_list))]
    if sum(1 for x in rr if x >= 2.0) >= obs_hi:
        cnt_hi += 1
    if max(rr) >= obs_max:
        cnt_max += 1
p_hi = float(cnt_hi) / (NPERM + 1)
p_max = float(cnt_max) / (NPERM + 1)

print()
print("=== NULL: IF RECOVERY TRACKED DATABASE SIZE ===")
print("  %d wild SGBs distributed across %d genera in proportion to their"
      % (tot, len(gen_list)))
print("  GTDB species-cluster counts, %d draws." % NPERM)
print("  observed genera with ratio >= 2.0 : %d   p = %.4f" % (obs_hi, p_hi))
print("  observed max ratio                : %.2f  p = %.4f" % (obs_max, p_max))
print()
print("  This null asks: if we had recovered SGBs in proportion to how well")
print("  each genus is already sequenced, would we still see expansions this")
print("  large? A small p means the expansion is not explained by that model.")

with open(OUTG, "w") as f:
    f.write("genus\twild_sgbs\tgtdb_species_clusters\tratio\n")
    for g, c, nr, x in grows:
        f.write("%s\t%d\t%d\t%s\n" % (g, c, nr, "%.3f" % x if nr else "NA"))
print()
print("wrote", OUTG)

print()
print("=== CALIBRATION AGAINST PUBLISHED CATALOGS ===")
print("  Youngblut 2020 mSystems : 1,522 SGBs, 1,184 novel species,")
print("                            266 novel genera, 6 novel families")
print("  Almeida 2021 Nat Biotech: 81%% of UHGG species had no cultured rep")
print("  Report YOUR percentages next to these. Novelty proportion is the")
print("  field-standard statistic and needs no arbitrary threshold.")
print()
print("REMINDER: the <=20 GTDB genomes / 87.2%% statistic is DEAD. Background")
print("is 88.3%% of genera, Mann-Whitney p = 0.470. Do not use it.")
print("DONE_NOVELTY_PROPORTIONS")
# SENTINEL_END
