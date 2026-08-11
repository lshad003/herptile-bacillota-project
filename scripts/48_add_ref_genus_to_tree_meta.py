# GTDB genus assignments are joined onto the reference tips, which carry no genus in the tree metadata.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/add_ref_genus_to_tree_meta.py
# Output: work/rep_tree/figure_tree_metadata_genus.tsv
# ADD_REF_GENUS_TO_TREE_META_V1_20260806
# figure_tree_metadata.tsv has a blank genus for all 1,247 gtdb_ref tips,
# so the tree cannot be collapsed to genus for 75% of its tips. This joins
# GTDB r220 taxonomy onto them and writes a NEW file. The original is not
# modified.
# Tip format is  gtdb_ref|gtdb_ref|GCA_000158655.1
# Taxonomy keys are RS_GCF_... and GB_GCA_..., so the prefix is stripped,
# the same way scripts/four_catalog_genus_table.py does it.

import os, sys
from collections import Counter

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
META = os.path.join(ROOT, "work/rep_tree/figure_tree_metadata.tsv")
OUT = os.path.join(ROOT, "work/rep_tree/figure_tree_metadata_genus.tsv")
GTDB_TAX = "/srv/projects/db/gtdbtk/220/taxonomy/gtdb_taxonomy.tsv"

if os.path.exists(OUT):
    raise SystemExit("REFUSING TO OVERWRITE %s, move it first" % OUT)
for p in (META, GTDB_TAX):
    if not os.path.exists(p):
        sys.exit("MISSING: %s" % p)


def parse_tax(s):
    out = {}
    for part in s.strip().split(";"):
        part = part.strip()
        if len(part) > 3 and part[1:3] == "__":
            out[part[0]] = part[3:]
    return out


ref_genus = {}
ref_family = {}
with open(GTDB_TAX) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        acc = f[0].strip().replace("GB_", "").replace("RS_", "")
        t = parse_tax(f[1])
        ref_genus[acc] = t.get("g", "")
        ref_family[acc] = t.get("f", "")
print("taxonomy accessions loaded: %d" % len(ref_genus))

with open(META) as fh:
    header = fh.readline().rstrip("\n").split("\t")
rows = []
with open(META) as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < len(header):
            f = f + [""] * (len(header) - len(f))
        rows.append(dict(zip(header, f)))
print("metadata rows: %d" % len(rows))

gi = header.index("genus")
fi = header.index("family")

filled = 0
unmatched = []
for r in rows:
    if r["arm"] != "gtdb_ref":
        continue
    if r["genus"].strip():
        continue
    acc = r["tip"].split("|")[-1].strip()
    g = ref_genus.get(acc)
    if g is None:
        unmatched.append(acc)
        continue
    r["genus"] = g
    if not r["family"].strip():
        r["family"] = ref_family.get(acc, "")
    filled += 1

print("gtdb_ref tips filled: %d" % filled)
print("gtdb_ref tips unmatched: %d" % len(unmatched))
for a in unmatched[:10]:
    print("   %s" % a)

blank = [r for r in rows if not r["genus"].strip()]
print()
print("tips still without a genus: %d" % len(blank))
print("  by arm: %s" % dict(Counter(r["arm"] for r in blank)))
print("  (herptile and other query tips with no GTDB genus call are expected)")

print()
print("distinct genera per arm:")
for arm in sorted(set(r["arm"] for r in rows)):
    gs = set(r["genus"] for r in rows if r["arm"] == arm and r["genus"].strip())
    n = sum(1 for r in rows if r["arm"] == arm)
    print("  %-16s %5d tips, %4d named genera" % (arm, n, len(gs)))

allg = set(r["genus"] for r in rows if r["genus"].strip() and r["arm"] != "outgroup")
print()
print("total distinct genera across non-outgroup tips: %d" % len(allg))

with open(OUT, "w") as f:
    f.write("\t".join(header) + "\n")
    for r in rows:
        f.write("\t".join(r.get(h, "") for h in header) + "\n")
print()
print("wrote %s" % OUT)
print("ADD_REF_GENUS_V1_20260806_COMPLETE")
