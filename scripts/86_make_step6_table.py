#!/usr/bin/env python3
# Supplementary table S10 is built for section 3.6.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/make_step6_table.py
# Output: tables/TableS10_bgc_class_composition.tsv
"""
Supplementary table S10 for section 3.6, product class composition by arm.
Refuses to overwrite an existing file.
"""
import os
import sys

AG = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
RP = "/bigdata/stajichlab/lshad003/herptile-bacillota-project"

SRC = os.path.join(AG, "results/bgc_class_composition.tsv")
OUT = os.path.join(RP, "tables/TableS10_bgc_class_composition.tsv")

ARMS = ["UHM amphibian", "EHI newt", "EHI mammal", "Youngblut", "GTDB reference"]

if os.path.exists(OUT):
    sys.exit("REFUSING TO OVERWRITE %s" % OUT)

rows = []
with open(SRC) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ix = {n: i for i, n in enumerate(hdr)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        rows.append({
            "cls": f[ix["product_class"]],
            "arm": f[ix["arm"]],
            "n": int(f[ix["n_genomes"]]),
            "prop": float(f[ix["mean_per_genome_proportion"]]),
            "flag": int(f[ix["n_arms_fragmentation_flagged"]]),
        })

print("rows: %d" % len(rows))
classes = sorted(set(r["cls"] for r in rows))
print("product classes: %d" % len(classes))

denom = {}
for r in rows:
    denom.setdefault(r["arm"], set()).add(r["n"])
print("")
print("DENOMINATOR PER ARM (genomes with >= 2 regions and an N50)")
for a in ARMS:
    d = denom.get(a, set())
    print("  %-18s %s%s" % (a, sorted(d),
                            "  CONSTANT" if len(d) == 1 else "  VARIES, CHECK"))

print("")
print("%-32s %6s %s" % ("product class", "flagged", "mean proportion by arm"))
byclass = {}
for r in rows:
    byclass.setdefault(r["cls"], {})[r["arm"]] = r
for c in sorted(classes, key=lambda x: byclass[x][ARMS[0]]["flag"]):
    d = byclass[c]
    vals = " ".join("%s %.3f" % (a.split()[-1][:4], d[a]["prop"])
                    for a in ARMS if a in d)
    print("%-32s %6d %s" % (c, d[ARMS[0]]["flag"], vals))

clean = sorted(c for c in classes if byclass[c][ARMS[0]]["flag"] == 0)
print("")
print("CLASSES PASSING THE GATE AT ZERO FLAGGED ARMS: %d" % len(clean))
for c in clean:
    d = byclass[c]
    lo = min((d[a]["prop"], a) for a in ARMS if a in d)
    hi = max((d[a]["prop"], a) for a in ARMS if a in d)
    print("  %-30s range %.3f (%s) to %.3f (%s)" % (c, lo[0], lo[1], hi[0], hi[1]))

with open(OUT, "w") as fh:
    fh.write("product_class\tarm\tn_genomes\tmean_per_genome_proportion\t"
             "n_arms_fragmentation_flagged\tpasses_gate\n")
    for c in sorted(classes):
        for a in ARMS:
            if a not in byclass[c]:
                continue
            r = byclass[c][a]
            fh.write("%s\t%s\t%d\t%.6f\t%d\t%s\n"
                     % (c, a, r["n"], r["prop"], r["flag"],
                        "yes" if r["flag"] == 0 else "no"))

print("")
print("WROTE %s" % OUT)
print("")
print("CAPTION REQUIREMENTS:")
print("  Composition is each genome's share of regions in a class, averaged")
print("  within arm, over Ruminococcaceae genomes carrying at least two")
print("  antiSMASH regions with an assembly N50 available. Classes were")
print("  included only where at least 50 regions were detected across all")
print("  arms, and an arm only where at least 20 genomes qualified.")
print("  The Youngblut arm sits at exactly 20 genomes, the inclusion floor,")
print("  so its means are the least stable in the table.")
print("  A class is flagged where its proportion tracks contiguity within an")
print("  arm, since a between-arm difference in such a class cannot be")
print("  separated from fragmentation.")
# MAKE_STEP6_TABLE_V1
