#!/usr/bin/env python3
# Supplementary table S6 is built for section 3.3.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/make_step3_table.py
# Output: tables/TableS6_four_catalog_genus_recovery.tsv
"""
Supplementary table S6 for section 3.3, genus recovery across four catalogs.
Refuses to overwrite an existing file.
"""
import os
import sys

AG = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
RP = "/bigdata/stajichlab/lshad003/herptile-bacillota-project"

SRC = os.path.join(AG, "results/four_catalog_genus_table.tsv")
OUT = os.path.join(RP, "tables/TableS6_four_catalog_genus_recovery.tsv")

if os.path.exists(OUT):
    sys.exit("REFUSING TO OVERWRITE %s" % OUT)

rows = []
with open(SRC) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ix = {n: i for i, n in enumerate(hdr)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        rows.append({
            "genus": f[ix["genus"]],
            "block": f[ix["block"]],
            "label": f[ix["block_label"]],
            "herptile": int(f[ix["herptile"]]),
            "ehi": int(f[ix["ehi"]]),
            "youngblut": int(f[ix["youngblut"]]),
            "gtdb": int(f[ix["gtdb_ref"]]),
        })

print("data rows: %d" % len(rows))
print("")
print("BLOCKS")
blocks = {}
for r in rows:
    blocks.setdefault(r["block"], []).append(r)
for b in sorted(blocks):
    s = blocks[b]
    print("  %-34s %3d genera, %4d herptile units, %s"
          % (b, len(s), sum(r["herptile"] for r in s), s[0]["label"]))

print("")
print("COLUMN TOTALS")
for col in ("herptile", "ehi", "youngblut", "gtdb"):
    print("  %-12s %d" % (col, sum(r[col] for r in rows)))
print("")
print("CHECK: the herptile column should total 193 if these are pooled dRep")
print("units and 220 if they are SGBs. The 3.3 text says 97 units in the")
print("15-genus block, so compare that against the block total above.")

rows.sort(key=lambda r: (r["block"], -r["herptile"], r["genus"]))
with open(OUT, "w") as fh:
    fh.write("genus\tblock\tblock_label\therptile_units\tehi_mammal_units\t"
             "youngblut_units\tgtdb_species_clusters\n")
    for r in rows:
        fh.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\n" % (
            r["genus"], r["block"], r["label"], r["herptile"],
            r["ehi"], r["youngblut"], r["gtdb"]))

print("")
print("WROTE %s, %d rows" % (OUT, len(rows)))
print("CAPTION: the ehi column is EHI mammal and bird only, not EHI newt.")
print("Genus-unassigned units are not represented in this table.")
# MAKE_STEP3_TABLE_V1
