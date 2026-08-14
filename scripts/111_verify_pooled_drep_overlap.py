#!/usr/bin/env python3
# Pooled dereplication clusters checked for any cluster containing a wild genome and a genome from another set.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/verify_pooled_drep_zero_overlap.py
# Output: printed report
"""
Verifies the Step 2 claim that no pooled 95% ANI cluster contains a wild
herptile genome together with a genome from any other set.

Reads Cdb.csv from the five-arm pooled dRep directly. Genome names carry
the staging prefixes herp (wild herptile), amph (EHI newt), ehi (EHI
mammal), yb (Youngblut), ref (GTDB reference).

Prints:
  - genomes and clusters by arm
  - number of clusters containing more than one arm (the 57 figure)
  - EVERY multi-arm cluster that contains a herptile genome, listed in
    full, or an explicit statement that there are none
  - the same for EHI newt, since the amphibian replication claim depends
    on the newt arm too

Read-only; prints only.
"""

import csv
import os
import sys
from collections import Counter, defaultdict

CDB = ("/bigdata/stajichlab/lshad003/ruminococcaceae-agent/"
       "work/pooled_drep_rum/drep_out/data_tables/Cdb.csv")

PREFIXES = {"herp": "wild herptile", "amph": "EHI newt",
            "ehi": "EHI mammal", "yb": "Youngblut", "ref": "GTDB reference"}

if not os.path.isfile(CDB):
    sys.exit("STOP: missing input:\n" + CDB)

members = defaultdict(list)
arm_counts = Counter()
unparsed = []

with open(CDB) as fh:
    rd = csv.DictReader(fh)
    if "secondary_cluster" not in rd.fieldnames:
        sys.exit("STOP: no secondary_cluster column; header is:\n"
                 + ",".join(rd.fieldnames))
    for r in rd:
        g = r["genome"].strip()
        if "__" not in g:
            unparsed.append(g)
            continue
        prefix = g.split("__", 1)[0]
        if prefix not in PREFIXES:
            unparsed.append(g)
            continue
        arm_counts[prefix] += 1
        members[r["secondary_cluster"].strip()].append((prefix, g))

print("genomes clustered: %d" % sum(arm_counts.values()))
for p in ("herp", "amph", "ehi", "yb", "ref"):
    print("  %-16s %d" % (PREFIXES[p], arm_counts[p]))
print("clusters: %d" % len(members))

if unparsed:
    print("")
    print("UNPARSED genome names: %d" % len(unparsed))
    for g in unparsed[:10]:
        print("  ", g)
    sys.exit("STOP: genome names outside the expected prefix scheme.")

multi = {}
for c, mem in members.items():
    arms = set(p for p, g in mem)
    if len(arms) > 1:
        multi[c] = (arms, mem)

print("")
print("clusters containing more than one set: %d" % len(multi))

combos = Counter()
for c, (arms, mem) in multi.items():
    combos[tuple(sorted(arms))] += 1
print("by combination:")
for combo, n in combos.most_common():
    print("  %3d  %s" % (n, " + ".join(PREFIXES[a] for a in combo)))

for focus, label in (("herp", "WILD HERPTILE"), ("amph", "EHI NEWT")):
    hits = {c: v for c, v in multi.items() if focus in v[0]}
    print("")
    print("=" * 62)
    print("%s genomes sharing a cluster with any other set: %d cluster(s)"
          % (label, len(hits)))
    if not hits:
        print("  NONE. The zero-overlap claim holds for this set.")
    else:
        for c, (arms, mem) in sorted(hits.items()):
            print("  cluster %s: %s" % (c, ", ".join(sorted(arms))))
            for p, g in sorted(mem):
                print("      %-16s %s" % (PREFIXES[p], g))

print("")
print("READ IT THIS WAY: the Step 2 sentence claims no cluster contains a")
print("wild-catalog genome and a genome from any other set. That claim is")
print("verified only if the WILD HERPTILE block above says NONE.")

# VERIFY_POOLED_DREP_ZERO_OVERLAP_V1
