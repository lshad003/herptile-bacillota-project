#!/usr/bin/env python3
# Supplementary tables S4 and S5 are built for section 3.2.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/make_step2_tables.py
# Output: tables/TableS4_genus_expansion.tsv, tables/TableS5_multiarm_clusters.tsv
"""
Supplementary tables for section 3.2.
S4 is per-genus expansion, S5 is the multi-arm cluster positive control.
Refuses to overwrite existing files.
"""
import os
import sys

AG = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
RP = "/bigdata/stajichlab/lshad003/herptile-bacillota-project"

EXP = os.path.join(AG, "results/genus_expansion.tsv")
MUL = os.path.join(AG, "results/pooled_drep_multiarm_audit.tsv")
S4 = os.path.join(RP, "tables/TableS4_genus_expansion.tsv")
S5 = os.path.join(RP, "tables/TableS5_multiarm_clusters.tsv")

for p in (S4, S5):
    if os.path.exists(p):
        sys.exit("REFUSING TO OVERWRITE %s" % p)

rows = []
with open(EXP) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ix = {n: i for i, n in enumerate(hdr)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        rows.append({
            "genus": f[ix["genus"]],
            "wild": int(f[ix["wild_sgbs"]]),
            "gtdb": int(f[ix["gtdb_species_clusters"]]),
            "ratio": float(f[ix["ratio"]]),
        })

rows.sort(key=lambda r: (-r["ratio"], r["genus"]))
total_wild = sum(r["wild"] for r in rows)
at2 = [r for r in rows if r["ratio"] >= 2.0]
exact2 = [r for r in rows if abs(r["ratio"] - 2.0) < 1e-9]
thin = [r for r in at2 if r["gtdb"] <= 2]

with open(S4, "w") as fh:
    fh.write("genus\twild_sgbs\tgtdb_r220_species_clusters\tratio\t"
             "at_least_twofold\tgtdb_clusters_le_2\n")
    for r in rows:
        fh.write("%s\t%d\t%d\t%.3f\t%s\t%s\n" % (
            r["genus"], r["wild"], r["gtdb"], r["ratio"],
            "yes" if r["ratio"] >= 2.0 else "no",
            "yes" if r["gtdb"] <= 2 else "no"))

print("TableS4: %d genera, %d wild SGBs" % (len(rows), total_wild))
print("  at or above twofold: %d" % len(at2))
print("  exactly twofold    : %s" % ", ".join(r["genus"] for r in exact2))
print("  resting on <=2 GTDB clusters: %d of %d" % (len(thin), len(at2)))
for t in (1.5, 2.0, 2.5, 3.0):
    print("  threshold %.1f -> %d genera"
          % (t, sum(1 for r in rows if r["ratio"] >= t)))

out = []
with open(MUL) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ix = {n: i for i, n in enumerate(hdr)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        out.append({
            "cluster": f[ix["secondary_cluster"]],
            "n_genomes": int(f[ix["n_genomes"]]),
            "n_arms": int(f[ix["n_arms"]]),
            "arms": f[ix["arms"]],
            "ids": f[ix["genome_ids"]],
        })

combos = {}
for r in out:
    combos[r["arms"]] = combos.get(r["arms"], 0) + 1
herp = [r for r in out if "herptile" in r["arms"] or "ehi_amphibian" in r["arms"]]

out.sort(key=lambda r: (-r["n_arms"], r["arms"], r["cluster"]))
with open(S5, "w") as fh:
    fh.write("secondary_cluster\tn_genomes\tn_arms\tarms\tgenome_ids\n")
    for r in out:
        fh.write("%s\t%d\t%d\t%s\t%s\n" % (
            r["cluster"], r["n_genomes"], r["n_arms"], r["arms"], r["ids"]))

print("")
print("TableS5: %d multi-arm clusters" % len(out))
for k in sorted(combos, key=lambda x: -combos[x]):
    print("  %-34s %d" % (k, combos[k]))
print("  clusters containing herptile or ehi_amphibian: %d" % len(herp))
if herp:
    print("  WARNING: THE ZERO CLAIM DOES NOT HOLD IN THIS FILE")

print("")
print("CAPTION REQUIREMENTS:")
print("  S4: ratio is wild SGBs divided by GTDB r220 species clusters of the")
print("      same genus. The twofold threshold was chosen after inspecting")
print("      the data. Genus-unassigned wild SGBs are not represented.")
print("  S5: these are the clusters demonstrating that joint dereplication")
print("      recovers cross-arm identity where it exists. Arm labels are")
print("      staged prefixes, not host categories.")
# MAKE_STEP2_TABLES_V1
