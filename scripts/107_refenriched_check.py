#!/usr/bin/env python3
# Reference-enriched orthologous groups annotated from member proteins.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/angelakisella_refenriched_check.py
# Output: printed report
"""
Sanity check on the reference-enriched side of the Angelakisella happi
result: annotates the top N reference-enriched OGs (diff < 0, ranked by
q) from actual member proteins in the non-amphibian genomes, so their
identity comes from real annotations, not OG-level description transfer.

Purpose: determine whether the near-universal genes at amphibian
prevalence zero (COG0504/pyrG, COG1490, COG0337, ...) are the known
pyrG-class amphibian-side absences (biologically real, already
verified in the project) or a detection artifact.

Also flags, per OG, whether ANY amphibian genome carries it at all
(prev_amphibian > 0), since a scattering of amphibian carriers argues
against a hard lineage absence.

Read-only; prints only.
"""

import os
import sys
from collections import Counter, defaultdict

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
HAPPI = B + "/results/happi_angelakisella_og.tsv"
M_G = B + "/results/angelakisella_matrix_genomes_v2.tsv"
EGG = B + "/results/eggnog"
EGG_REFS = B + "/results/eggnog_refs"

TOP_N = 20
N_GENOMES_SCAN = 6   # reference genomes to pull annotations from

for p in (HAPPI, M_G):
    if not os.path.isfile(p):
        sys.exit("STOP: missing input:\n" + p)

rows = []
with open(HAPPI) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        f = dict(zip(hdr, line.rstrip("\n").split("\t")))
        if float(f["qvalue"]) < 0.05 and float(f["diff"]) < 0:
            rows.append(f)

rows.sort(key=lambda r: float(r["pvalue"]))
top = rows[:TOP_N]
print("reference-enriched significant OGs: %d, annotating top %d"
      % (len(rows), len(top)))

nonamph = []
with open(M_G) as fh:
    fh.readline()
    for line in fh:
        g, cat, grp = line.rstrip("\n").split("\t")
        if grp == "non_amphibian":
            nonamph.append(g)

scan_genomes = nonamph[:N_GENOMES_SCAN]
want = {r["og"] for r in top}

ann = defaultdict(lambda: {"ko": Counter(), "name": Counter(),
                           "desc": Counter()})
for key in scan_genomes:
    arm, gid = key.split("__", 1)
    d = EGG_REFS if arm == "gtdb_ref" else EGG
    path = os.path.join(d, gid + ".emapper.annotations")
    if not os.path.isfile(path):
        sys.exit("STOP: missing annotations: " + path)
    hdr2 = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#query"):
                hdr2 = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if hdr2 is None or len(f) < len(hdr2):
                continue
            row = dict(zip(hdr2, f))
            for og in row.get("eggNOG_OGs", "").split(","):
                og = og.strip()
                if "@" not in og:
                    continue
                unit, tax = og.split("@", 1)
                if tax.split("|", 1)[0] == "2" and unit in want:
                    ann[unit]["ko"][row.get("KEGG_ko", "-")] += 1
                    ann[unit]["name"][row.get("Preferred_name", "-")] += 1
                    ann[unit]["desc"][row.get("Description", "-")[:70]] += 1

print("annotations pulled from %d reference genomes: %s"
      % (len(scan_genomes), ", ".join(scan_genomes)))
print("")
print("%-9s %6s %6s %8s  %s" % ("OG", "prevA", "prevN", "q", "identity"))
for r in top:
    og = r["og"]
    a = ann.get(og)
    if a and a["name"]:
        name = a["name"].most_common(1)[0][0]
        ko = a["ko"].most_common(1)[0][0]
        desc = a["desc"].most_common(1)[0][0]
        ident = "%s | %s | %s" % (name, ko, desc)
    else:
        ident = "(no member found in scanned genomes)"
    marker = ""
    if float(r["prev_amphibian"]) == 0.0:
        marker = "  <-- amphibian ZERO"
    print("%-9s %6.2f %6.2f %8.1e  %s%s"
          % (og, float(r["prev_amphibian"]),
             float(r["prev_non_amphibian"]),
             float(r["qvalue"]), ident, marker))

print("")
print("READ IT THIS WAY:")
print("Amphibian-ZERO rows that are pyrimidine synthesis (pyrG/COG0504)")
print("or its co-lost neighbors corroborate the project's verified")
print("amphibian-side absences and validate the happi run. Amphibian-ZERO")
print("rows that are random unrelated housekeeping genes would instead")
print("suggest a systematic amphibian-side detection problem and would")
print("need the same absence verification (HMM/tblastn) before the 866")
print("figure is quoted.")

# ANGELAKISELLA_REFENRICHED_CHECK_V1
