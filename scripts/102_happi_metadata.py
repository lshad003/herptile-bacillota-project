#!/usr/bin/env python3
# CheckM2 completeness joined onto the 57 test genomes.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/angelakisella_happi_metadata.py
# Output: results/angelakisella_happi_metadata.tsv
"""
Builds the metadata table for happi on the Angelakisella OG matrix:
the 57 Angelakisella genomes (31 amphibian, 26 non_amphibian; the 4
Heteroruminococcus context genomes are excluded from the test), with
CheckM2 completeness joined from the verified five-arm run
(work/checkm2_all_arms/checkm2_out_v2/quality_report.tsv).

CheckM2 Name format: prefix__ID.fa with prefixes
  herp=herptile, amph=ehi_amphibian, ref=gtdb_ref, yb=youngblut.
Split on the FIRST "__" only (Youngblut IDs contain "__").

Refuses to run on any missing or ambiguous completeness match.
Writes results/angelakisella_happi_metadata.tsv (refuses overwrite).
"""

import os
import sys

BASE = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
M_G = os.path.join(BASE, "results/angelakisella_matrix_genomes_v2.tsv")
QR = os.path.join(BASE,
                  "work/checkm2_all_arms/checkm2_out_v2/quality_report.tsv")
OUT = os.path.join(BASE, "results/angelakisella_happi_metadata.tsv")

ARM2PREFIX = {
    "herptile": "herp",
    "ehi_amphibian": "amph",
    "gtdb_ref": "ref",
    "youngblut": "yb",
}

if os.path.exists(OUT):
    sys.exit("STOP: output exists, refusing to overwrite:\n" + OUT)
for p in (M_G, QR):
    if not os.path.isfile(p):
        sys.exit("STOP: missing input:\n" + p)

comp = {}
with open(QR) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    i_name = header.index("Name")
    i_comp = header.index("Completeness")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        name = f[i_name]
        if "__" not in name:
            continue
        prefix, rest = name.split("__", 1)
        for ext in (".fa.gz", ".fna", ".fa"):
            if rest.endswith(ext):
                rest = rest[:-len(ext)]
                break
        key = (prefix, rest)
        if key in comp:
            sys.exit("STOP: duplicate CheckM2 row for %s__%s" % key)
        comp[key] = float(f[i_comp])

print("CheckM2 rows indexed:", len(comp))
if len(comp) != 1957:
    sys.exit("STOP: expected 1957 CheckM2 rows, got %d." % len(comp))

rows = []
with open(M_G) as fh:
    fh.readline()
    for line in fh:
        g, cat, grp = line.rstrip("\n").split("\t")
        if grp == "context":
            continue
        arm, gid = g.split("__", 1)
        prefix = ARM2PREFIX.get(arm)
        if prefix is None:
            sys.exit("STOP: unmapped catalog: " + arm)
        key = (prefix, gid)
        if key not in comp:
            sys.exit("STOP: no CheckM2 completeness for %s (looked up %s__%s)"
                     % (g, prefix, gid))
        rows.append((g, grp, comp[key]))

n_a = sum(1 for g, gr, c in rows if gr == "amphibian")
n_n = sum(1 for g, gr, c in rows if gr == "non_amphibian")
print("genomes: %d amphibian, %d non_amphibian" % (n_a, n_n))
if (n_a, n_n) != (31, 26):
    sys.exit("STOP: expected 31 + 26.")

with open(OUT, "w") as out:
    out.write("genome\tgroup\tcompleteness\n")
    for g, gr, c in sorted(rows):
        out.write("%s\t%s\t%.2f\n" % (g, gr, c))

ca = [c for g, gr, c in rows if gr == "amphibian"]
cn = [c for g, gr, c in rows if gr == "non_amphibian"]
print("completeness mean: amphibian %.1f, non_amphibian %.1f"
      % (sum(ca) / len(ca), sum(cn) / len(cn)))
print("wrote:", OUT)

# ANGELAKISELLA_HAPPI_METADATA_V1
