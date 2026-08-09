#!/usr/bin/env python3
# Summary of chimerism screening across the SGB representatives.
#
# Source: ruminococcaceae-agent/scripts/gunc_audit.py
# Reads:  results/gunc_sgb/GUNC.progenomes_2.1.maxCSS_level.tsv
#         data/sgb_manifest.tsv
# Writes: results/gunc_audit_by_sgb.tsv
#
# Two thresholds are applied. The first is GUNC's own pass criterion at a
# clade separation score of 0.45. The second is a study-defined criterion
# requiring a clade separation score of at least 0.85, a reference
# representation score of at least 0.5, and a maximum contamination signal at
# family level or above.
#
# Reference representation scores are low across this catalog, with a median
# of 0.48, so detection is conservative. Results are reported as no chimerism
# detected under limited reference representation rather than as no chimerism
# present.
import os, sys
from collections import Counter, defaultdict
import numpy as np

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
GUNC = B + "/results/gunc_sgb/GUNC.progenomes_2.1.maxCSS_level.tsv"
SGB = B + "/data/sgb_manifest.tsv"
OUT = B + "/results/gunc_audit_by_sgb.tsv"

FAM = "Ruminococcaceae"
FOCAL = ("UBA866", "Anaerotruncus", "Angelakisella", "Ruthenibacterium")
RRS_MIN = 0.5
CSS_HARD = 0.85

for p in (GUNC, SGB):
    if not os.path.exists(p):
        print("MISSING:", p); sys.exit(1)

g = {}
with open(GUNC) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    I = {k: i for i, k in enumerate(h)}
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(I.values()):
            continue
        try:
            g[p[I["genome"]]] = dict(
                lvl=p[I["taxonomic_level"]],
                css=float(p[I["clade_separation_score"]]),
                cont=float(p[I["contamination_portion"]]),
                rrs=float(p[I["reference_representation_score"]]),
                ident=float(p[I["mean_hit_identity"]]),
                ngenes=int(p[I["n_genes_called"]]),
                nmap=int(p[I["n_genes_mapped"]]),
                ncontig=int(p[I["n_contigs"]]),
                surplus=float(p[I["n_effective_surplus_clades"]]),
                passed=p[I["pass.GUNC"]].strip().lower() == "true")
        except ValueError:
            continue
print("GUNC rows: %d" % len(g))

meta = {}
with open(SGB) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    J = {k: i for i, k in enumerate(h)}
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(J.values()):
            continue
        meta[p[J["representative"]]] = dict(
            sgb=p[J["sgb"]], fam=p[J["family"]], gen=p[J["genus"]],
            wild=p[J["has_wild"]] == "yes",
            comp=float(p[J["rep_completeness"]]),
            cont=float(p[J["rep_contamination"]]),
            n=int(p[J["n_mags"]]))
print("SGB manifest: %d" % len(meta))

both = [k for k in g if k in meta]
print("joined: %d" % len(both))
if len(both) < 0.95 * len(g):
    print("JOIN PROBLEM"); sys.exit(1)

print()
print("=== DEFAULT GUNC VERDICT, ALL %d SGBs ===" % len(both))
npass = sum(1 for k in both if g[k]["passed"])
print("  pass.GUNC True  : %d (%.1f%%)" % (npass, 100.0 * npass / len(both)))
print("  pass.GUNC False : %d (%.1f%%)"
      % (len(both) - npass, 100.0 * (len(both) - npass) / len(both)))
print("  (published human gut MAG sets: 15-30%% flagged after CheckM filtering)")

print()
print("=== TAXONOMIC LEVEL OF MAX CSS ===")
for lvl, n in Counter(g[k]["lvl"] for k in both).most_common():
    fl = sum(1 for k in both if g[k]["lvl"] == lvl and not g[k]["passed"])
    print("  %-12s %5d SGBs | %4d fail" % (lvl, n, fl))

rrs = np.array([g[k]["rrs"] for k in both])
print()
print("=== REFERENCE REPRESENTATION ===")
print("  RRS  min %.3f  median %.3f  max %.3f" % (rrs.min(), np.median(rrs), rrs.max()))
print("  RRS < %.1f (poorly represented) : %d (%.1f%%)"
      % (RRS_MIN, int((rrs < RRS_MIN).sum()), 100.0 * (rrs < RRS_MIN).sum() / len(rrs)))
low = [k for k in both if g[k]["rrs"] < RRS_MIN]
print("  of those, pass.GUNC False       : %d"
      % sum(1 for k in low if not g[k]["passed"]))
print("  A low RRS means the lineage is undersampled in proGenomes, so a bad")
print("  score there is weak evidence of real chimerism.")

def strict(k):
    return (g[k]["css"] >= CSS_HARD and g[k]["rrs"] >= RRS_MIN
            and g[k]["lvl"] in ("kingdom", "phylum", "class", "order", "family"))

st = [k for k in both if strict(k)]
print()
print("=== STRICT CALL (CSS >= %.2f, RRS >= %.1f, family level or above) ===" % (CSS_HARD, RRS_MIN))
print("  putative chimeras: %d (%.1f%%)" % (len(st), 100.0 * len(st) / len(both)))
print("  This is the conservative standard used where reference coverage is thin.")

def rep(sel, lab):
    if not sel:
        print("\n  %s: none" % lab); return
    print()
    print("  --- %s (%d) ---" % (lab, len(sel)))
    print("  %-28s %-8s %6s %6s %6s %6s %5s %5s"
          % ("genome", "level", "CSS", "cont", "RRS", "ident", "comp", "nMAG"))
    for k in sorted(sel, key=lambda z: -g[z]["css"])[:20]:
        m = meta[k]
        print("  %-28s %-8s %6.2f %6.2f %6.2f %6.2f %5.1f %5d"
              % (k[:28], g[k]["lvl"][:8], g[k]["css"], g[k]["cont"],
                 g[k]["rrs"], g[k]["ident"], m["comp"], m["n"]))

rep([k for k in st if meta[k]["fam"] == FAM], "%s, strict chimera calls" % FAM)

print()
print("=== %s SGBs ===" % FAM)
rs = [k for k in both if meta[k]["fam"] == FAM]
rw = [k for k in rs if meta[k]["wild"]]
print("  total %d | wild %d" % (len(rs), len(rw)))
print("  pass.GUNC False        : %d (%.1f%%)"
      % (sum(1 for k in rs if not g[k]["passed"]),
         100.0 * sum(1 for k in rs if not g[k]["passed"]) / len(rs)))
print("  strict chimera calls   : %d" % sum(1 for k in rs if strict(k)))
print("  RRS < %.1f              : %d" % (RRS_MIN, sum(1 for k in rs if g[k]["rrs"] < RRS_MIN)))

print()
print("=== FOCAL GENERA ===")
print("  %-18s %5s %8s %8s %8s %8s"
      % ("genus", "SGBs", "failGUNC", "strict", "RRS<0.5", "medCSS"))
for gg in FOCAL:
    sel = [k for k in rs if meta[k]["gen"] == gg]
    if not sel:
        continue
    css = sorted(g[k]["css"] for k in sel)
    print("  %-18s %5d %8d %8d %8d %8.3f"
          % (gg, len(sel),
             sum(1 for k in sel if not g[k]["passed"]),
             sum(1 for k in sel if strict(k)),
             sum(1 for k in sel if g[k]["rrs"] < RRS_MIN),
             css[len(css) // 2]))

fl = [k for k in both if not g[k]["passed"]]
if fl:
    cf = np.array([meta[k]["comp"] for k in fl])
    cp = np.array([meta[k]["comp"] for k in both if g[k]["passed"]])
    xf = np.array([meta[k]["cont"] for k in fl])
    xp = np.array([meta[k]["cont"] for k in both if g[k]["passed"]])
    print()
    print("=== DOES CheckM SEE WHAT GUNC SEES? ===")
    print("  completeness  fail %.2f | pass %.2f" % (cf.mean(), cp.mean()))
    print("  contamination fail %.2f | pass %.2f" % (xf.mean(), xp.mean()))
    print("  CheckM contamination >5%% among GUNC failures: %d of %d"
          % (int((xf > 5).sum()), len(fl)))
    print("  GUNC detects chimerism CheckM misses, so overlap is expected to be low.")

with open(OUT, "w") as f:
    f.write("genome\tsgb\tfamily\tgenus\thas_wild\tn_mags\tcheckm_completeness\t"
            "checkm_contamination\ttaxonomic_level\tclade_separation_score\t"
            "contamination_portion\treference_representation_score\t"
            "mean_hit_identity\tn_genes_called\tn_genes_mapped\tn_contigs\t"
            "pass_gunc\tstrict_chimera\n")
    for k in sorted(both):
        m = meta[k]; v = g[k]
        f.write("%s\t%s\t%s\t%s\t%s\t%d\t%.2f\t%.2f\t%s\t%.3f\t%.3f\t%.3f\t"
                "%.3f\t%d\t%d\t%d\t%s\t%s\n"
                % (k, m["sgb"], m["fam"], m["gen"], "yes" if m["wild"] else "no",
                   m["n"], m["comp"], m["cont"], v["lvl"], v["css"], v["cont"],
                   v["rrs"], v["ident"], v["ngenes"], v["nmap"], v["ncontig"],
                   "pass" if v["passed"] else "FAIL",
                   "yes" if strict(k) else "no"))
print()
print("wrote", OUT)
print()
print("TWO NUMBERS FOR THE PAPER: the default pass.GUNC rate, and the strict")
print("call requiring RRS >= 0.5 and family level or above. Report both. The")
print("strict set is what should actually be excluded.")
print("DONE_GUNC_AUDIT")
# SENTINEL_END
