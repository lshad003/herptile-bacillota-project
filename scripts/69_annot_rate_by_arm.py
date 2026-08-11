#!/usr/bin/env python3
# Per-protein annotation rate is compared between arms on both annotation layers.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/annot_rate_by_arm_v2.py
# Output: results/annot_rate_by_arm_v2.tsv
import os, sys
from scipy.stats import mannwhitneyu, spearmanr

BASE = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
SUM  = os.path.join(BASE, "results/annot_154_threshold_summary.tsv")
QUAL = os.path.join(BASE, "work/focal_genus_pangenome/checkm2_out/quality_report.tsv")
OUT  = os.path.join(BASE, "results/annot_rate_by_arm_v2.tsv")

def die(m):
    sys.stderr.write("FATAL: " + m + "\n"); sys.exit(1)

if os.path.exists(OUT):
    die("output exists, refusing to overwrite: " + OUT)
for p in (SUM, QUAL):
    if not os.path.isfile(p):
        die("missing " + p)

with open(SUM) as fh:
    head = fh.readline().rstrip("\n").split("\t")
    recs = [dict(zip(head, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]
for c in ("genome_id", "set", "arm", "genus", "kofam_annot_rate", "pfam_annot_rate", "n_proteins"):
    if c not in head:
        die("column %s absent from %s" % (c, SUM))
if len(recs) != 154:
    die("expected 154 summary rows, found %d" % len(recs))

# CheckM2 Name is Genus__arm__genomeid. Take the field after the LAST '__'.
qual = {}
with open(QUAL) as fh:
    qh = fh.readline().rstrip("\n").split("\t")
    for c in ("Name", "Completeness", "Contig_N50", "Total_Contigs"):
        if c not in qh:
            die("quality_report.tsv lacks column " + c)
    ni = qh.index("Name"); ci = qh.index("Completeness")
    n5 = qh.index("Contig_N50"); tc = qh.index("Total_Contigs")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if len(f) <= max(ni, ci, n5, tc):
            continue
        name = f[ni]
        if "__" not in name:
            die("unexpected CheckM2 Name with no '__': " + name)
        gid = name.rsplit("__", 1)[1]
        if gid in qual:
            die("duplicate genome id after Name transform: " + gid)
        qual[gid] = (float(f[ci]), float(f[n5]), float(f[tc]))
if len(qual) != 125:
    die("Name transform yielded %d unique ids, expected 125" % len(qual))

m = [r for r in recs if r["set"] == "focal125" and r["arm"] in ("amphibian", "reference")]
if len(m) != 124:
    die("R10 filter gave %d genomes, expected 124" % len(m))

for r in m:
    g = r["genome_id"]
    if g not in qual:
        die("no CheckM2 row for " + g)
    r["comp"], r["n50"], r["ctg"] = qual[g]
    r["kofam"] = float(r["kofam_annot_rate"])
    r["pfam"]  = float(r["pfam_annot_rate"])
    r["npro"]  = int(r["n_proteins"])

amp = [r for r in m if r["arm"] == "amphibian"]
ref = [r for r in m if r["arm"] == "reference"]
print("design: %d amphibian, %d reference" % (len(amp), len(ref)))
print("")

def med(v):
    v = sorted(v); n = len(v)
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2.0

lines = []
print("=== ARM CONTRAST, all 124 ===")
for lab, key in (("kofam", "kofam"), ("pfam", "pfam"), ("completeness", "comp"),
                 ("contig_N50", "n50"), ("n_contigs", "ctg"), ("n_proteins", "npro")):
    a = [r[key] for r in amp]; b = [r[key] for r in ref]
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    print("  %-12s amphibian median %.4f | reference median %.4f | MWU p %.4g"
          % (lab, med(a), med(b), p))
    lines.append(["all", "arm_contrast", lab, len(a), len(b), "%.5f" % med(a), "%.5f" % med(b), "%.6g" % p])
print("")

print("=== ARM CONTRAST WITHIN GENUS ===")
for gen in sorted(set(r["genus"] for r in m)):
    ga = [r for r in amp if r["genus"] == gen]; gb = [r for r in ref if r["genus"] == gen]
    print("  %s: %d amphibian, %d reference" % (gen, len(ga), len(gb)))
    if len(ga) < 3 or len(gb) < 3:
        print("    skipped, fewer than 3 per side")
        continue
    for lab, key in (("kofam", "kofam"), ("pfam", "pfam"), ("completeness", "comp"), ("contig_N50", "n50")):
        a = [r[key] for r in ga]; b = [r[key] for r in gb]
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        print("    %-12s amph %.4f | ref %.4f | MWU p %.4g" % (lab, med(a), med(b), p))
        lines.append([gen, "within_genus", lab, len(a), len(b), "%.5f" % med(a), "%.5f" % med(b), "%.6g" % p])
print("")

print("=== ANNOTATION RATE vs QUALITY, which mechanism ===")
for lab, key in (("kofam", "kofam"), ("pfam", "pfam")):
    for qlab, qkey in (("completeness", "comp"), ("contig_N50", "n50")):
        for grp, rs in (("all124", m), ("amphibian", amp), ("reference", ref)):
            rho, p = spearmanr([r[key] for r in rs], [r[qkey] for r in rs])
            print("  %-6s vs %-12s %-10s rho %+.3f  p %.4g  n %d" % (lab, qlab, grp, rho, p, len(rs)))
            lines.append([grp, "spearman_" + qkey, lab, len(rs), "", "%.4f" % rho, "", "%.6g" % p])
print("")

with open(OUT, "w") as out:
    out.write("group\ttest\tlayer\tn_a\tn_b\tstat_a\tstat_b\tp\n")
    for r in lines:
        out.write("\t".join(str(x) for x in r) + "\n")
print("WROTE: " + OUT)
# ANNOT_RATE_BY_ARM_V2_20260808
