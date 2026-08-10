#!/usr/bin/env python3
# Genus assignments for the focal genera are mapped from GTDB r220 onto r226 taxonomy.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/uba866_r226_map.py
# Output: results/genus_r220_to_r226_map.tsv
import os, sys
from collections import Counter

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
T220 = "/srv/projects/db/gtdbtk/220/taxonomy/bac120_taxonomy_r220_reps.tsv"
T226 = "/srv/projects/db/gtdbtk/226/taxonomy/bac120_taxonomy_r226_reps.tsv"
OUT = B + "/results/genus_r220_to_r226_map.tsv"

WATCH = ["UBA866", "Muricomes", "Anaerotruncus", "Angelakisella",
         "Ruthenibacterium", "Fimivivens", "JAAYCI01", "Gemmiger",
         "Faecousia", "Limivicinus", "Otoolea"]

for p in (T220, T226):
    if not os.path.exists(p):
        print("MISSING:", p); sys.exit(1)

def fld(t, pre):
    for f in t.split(";"):
        f = f.strip()
        if f.startswith(pre):
            return f[len(pre):]
    return ""

def load(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                d[p[0]] = p[1]
    return d

a = load(T220); b = load(T226)
print("r220 accessions: %d | r226: %d" % (len(a), len(b)))
print("shared accessions: %d" % len(set(a) & set(b)))

print()
print("=== WHERE EACH r220 GENUS WENT IN r226 ===")
rows = []
for g in WATCH:
    accs = [k for k, v in a.items() if fld(v, "g__") == g]
    if not accs:
        print("%-18s not found in r220" % g)
        continue
    dest = Counter()
    gone = 0
    for k in accs:
        if k in b:
            dest[fld(b[k], "g__") or "(none)"] += 1
        else:
            gone += 1
    n226 = sum(1 for v in b.values() if fld(v, "g__") == g)
    print("%-18s r220 n=%-4d | genus name still in r226: %d genomes"
          % (g, len(accs), n226))
    if gone:
        print("%-18s   %d r220 genomes are NOT in r226 at all" % ("", gone))
    for d, c in dest.most_common():
        fam = ""
        for k in accs:
            if k in b and (fld(b[k], "g__") or "(none)") == d:
                fam = fld(b[k], "f__"); break
        flag = "" if d == g else "   <== RENAMED"
        print("%-18s   -> %-22s %3d  (f__%s)%s" % ("", d, c, fam, flag))
        rows.append((g, len(accs), n226, d, c, fam, gone))
    print()

print("=== SIZE OF THE DESTINATION GENUS IN r226 ===")
seen = set()
for r in rows:
    d = r[3]
    if d in seen or d == "(none)":
        continue
    seen.add(d)
    n = sum(1 for v in b.values() if fld(v, "g__") == d)
    print("  %-24s %d genomes in r226" % (d, n))

with open(OUT, "w") as f:
    f.write("r220_genus\tn_r220\tn_r226_same_name\tr226_genus\tn_moved\t"
            "r226_family\tn_dropped_from_r226\n")
    for r in rows:
        f.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\n"
                % (r[0], r[1], r[2], r[3], r[4], r[5], r[6]))
print()
print("wrote", OUT)
print()
print("A renamed genus means the r220 and r226 counts in candidate_genus_table")
print("are NOT comparable for that row. The denominator for the representation")
print("claim is the r226 destination genus, not the r220 name.")
print("DONE_R226_MAP")
# SENTINEL_END
