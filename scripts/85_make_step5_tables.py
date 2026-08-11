#!/usr/bin/env python3
# Supplementary tables S8 and S9 are built for section 3.5.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/make_step5_tables.py
# Output: tables/TableS8_genus_interleaving.tsv,
#         tables/TableS9_cazy_significant_families.tsv
"""
Supplementary tables S8 and S9 for section 3.5.
S8 is the genus interleaving test, S9 the significant CAZy families.
Refuses to overwrite existing files.
"""
import os
import sys

AG = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
RP = "/bigdata/stajichlab/lshad003/herptile-bacillota-project"

INT = os.path.join(AG, "results/genus_interleaving.tsv")
CAZ = os.path.join(AG, "results/happi_cazy_focal.tsv")
S8 = os.path.join(RP, "tables/TableS8_genus_interleaving.tsv")
S9 = os.path.join(RP, "tables/TableS9_cazy_significant_families.tsv")

MIN_CHANGES = 5
Q = 0.05

for p in (S8, S9):
    if os.path.exists(p):
        sys.exit("REFUSING TO OVERWRITE %s" % p)

rows = []
whole = None
with open(INT) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ix = {n: i for i, n in enumerate(hdr)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        rec = {k: (f[i] if i < len(f) else "") for k, i in ix.items()}
        if rec["genus"] == "WHOLE_TREE":
            whole = rec
            continue
        rows.append(rec)

print("genera tested for interleaving: %d" % len(rows))
if whole:
    print("WHOLE_TREE row held out: %s transitions, ratio %s, p %s"
          % (whole["obs_changes"], whole["ratio"], whole["p"]))
    print("  it is a tree-wide summary, not a genus, and is excluded from S8")

rows.sort(key=lambda r: -int(r["obs_changes"]))
sel = [r for r in rows if int(r["obs_changes"]) >= MIN_CHANGES]
print("")
print("genera at >= %d transitions: %d  (%s)"
      % (MIN_CHANGES, len(sel), ", ".join(r["genus"] for r in sel)))
mismatch = [r for r in rows
            if ("testable" in r["verdict"] and "NOT" not in r["verdict"])
            != (int(r["obs_changes"]) >= MIN_CHANGES)]
print("genera where the file verdict disagrees with the applied criterion: %d"
      % len(mismatch))
for r in mismatch:
    print("  %-20s verdict '%s', %s transitions, p %s"
          % (r["genus"], r["verdict"], r["obs_changes"], r["p"]))

with open(S8, "w") as fh:
    fh.write("genus\tn_amphibian\tn_other\tobserved_transitions\tnull_mean\t"
             "null_lo\tnull_hi\tobs_over_null\tp\tlargest_pure_amphibian_clade\t"
             "verdict_in_source\tused_for_gene_content_testing\n")
    for r in rows:
        used = "yes" if int(r["obs_changes"]) >= MIN_CHANGES else "no"
        fh.write("\t".join([
            r["genus"], r["n_amphibian"], r["n_other"], r["obs_changes"],
            r["null_mean"], r["null_lo"], r["null_hi"], r["ratio"], r["p"],
            r["largest_pure_amph_clade"], r["verdict"], used]) + "\n")

print("")
print("WROTE %s, %d rows" % (S8, len(rows)))

fams = []
with open(CAZ) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    cx = {n: i for i, n in enumerate(hdr)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[cx["q"]] in ("", "NA"):
            continue
        if float(f[cx["q"]]) >= Q:
            continue
        fams.append({
            "family": f[cx["family"]],
            "na": f[cx["n_amphibian"]], "nr": f[cx["n_reference"]],
            "pa": float(f[cx["prev_amphibian"]]),
            "pr": float(f[cx["prev_reference"]]),
            "beta": f[cx["beta"]], "lrt": f[cx["LRT"]],
            "p": f[cx["p"]], "q": f[cx["q"]],
        })

hi_a = [x for x in fams if x["pa"] > x["pr"]]
hi_r = [x for x in fams if x["pr"] > x["pa"]]
print("")
print("CAZy families at q < %.2f: %d  (amphibian-higher %d, reference-higher %d)"
      % (Q, len(fams), len(hi_a), len(hi_r)))

ratios = []
for x in hi_r:
    if x["pa"] > 0:
        ratios.append((x["pr"] / x["pa"], x["family"]))
if ratios:
    ratios.sort()
    print("smallest reference-higher prevalence ratio: %.2f-fold (%s)"
          % ratios[0])
    zero = [x["family"] for x in hi_r if x["pa"] == 0]
    if zero:
        print("reference-higher families absent from amphibian genomes: %s"
              % ", ".join(zero))

fams.sort(key=lambda x: (-(x["pa"] - x["pr"]), x["family"]))
with open(S9, "w") as fh:
    fh.write("family\tdirection\tn_amphibian_present\tn_reference_present\t"
             "prevalence_amphibian\tprevalence_reference\tbeta\tLRT\tp\tq\n")
    for x in fams:
        d = "amphibian_higher" if x["pa"] > x["pr"] else "reference_higher"
        fh.write("%s\t%s\t%s\t%s\t%.4f\t%.4f\t%s\t%s\t%s\t%s\n"
                 % (x["family"], d, x["na"], x["nr"], x["pa"], x["pr"],
                    x["beta"], x["lrt"], x["p"], x["q"]))

print("WROTE %s, %d rows" % (S9, len(fams)))
print("")
print("CAPTION REQUIREMENTS:")
print("  S8: the applied criterion is the observed transition count, not the")
print("      verdict string, which calls underpowered genera testable. The")
print("      source verdict is retained so the two can be compared.")
print("  S9: prevalences are proportions of genomes in each arm; the model is")
print("      presence ~ amphibian + genus with completeness as the quality")
print("      variable, so genus is adjusted out and the tested term is the")
print("      within-genus amphibian contrast.")
# MAKE_STEP5_TABLES_V1
