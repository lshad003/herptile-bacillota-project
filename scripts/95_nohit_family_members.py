#!/usr/bin/env python3
# Full membership, lengths and Prodigal partial flags of the single cluster without reference homologs.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/angelakisella_nohit_cluster_check.py
# Output: printed report
"""
Dissects the single amphibian-specific Angelakisella cluster whose
amphibian members had no DIAMOND hit to the 25 reference proteomes:
  rep = gtdb_ref__GCA_937894495.1|CALCRR010000048.1_3

Lists every member with arm, genome, length, and Prodigal partial flag,
checks whether GCA_937894495.1 is one of the 25 allowed references,
and reports any raw DIAMOND hits for these queries that the qcov filter
removed.

Read-only; prints only.
"""

import csv
import glob
import os
import re
import sys
from collections import defaultdict

BASE = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"

CLU = os.path.join(BASE, "work/all_arms_pangenome/mmseqs/all_arms_clusters.tsv")
PROTEOMES = os.path.join(BASE, "work/all_arms_pangenome/proteomes")
META = os.path.join(BASE, "work/rep_tree/figure_tree_metadata_genus.tsv")
HITS = os.path.join(BASE, "work/angelakisella_diamond/hits.tsv")

REP = "gtdb_ref__GCA_937894495.1|CALCRR010000048.1_3"

GENUS = "Angelakisella"
KEEP_ARMS = {"herptile", "ehi_amphibian", "gtdb_ref"}
AMPH_ARMS = {"herptile", "ehi_amphibian"}


def core_genome_id(x):
    x = x.strip()
    if "__" in x:
        prefix, rest = x.split("__", 1)
        if prefix in KEEP_ARMS:
            x = rest
    elif "|" in x:
        prefix, rest = x.split("|", 1)
        if prefix in KEEP_ARMS:
            x = rest
    if "|" in x:
        x = x.split("|", 1)[0]
    if x.startswith(("RS_", "GB_")):
        x = x[3:]
    return x


allowed = set()
with open(META) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["genus"].strip() != GENUS:
            continue
        arm = r["arm"].strip()
        if arm not in KEEP_ARMS:
            continue
        allowed.add((arm, core_genome_id(r["genome"])))

rep_gid = core_genome_id(REP.split("|", 1)[0].split("__", 1)[1])
print("rep genome:", rep_gid)
print("rep genome in 25 allowed references:",
      ("gtdb_ref", rep_gid) in allowed)


def parse_member(pid):
    if "|" not in pid:
        return None
    genome_prefix, local_pid = pid.split("|", 1)
    if "__" not in genome_prefix:
        return None
    arm, raw_gid = genome_prefix.split("__", 1)
    return arm, raw_gid, genome_prefix, local_pid


members = []
with open(CLU) as fh:
    for line in fh:
        rep, mem = line.rstrip("\n").split("\t")
        if rep == REP:
            members.append(mem)

print("")
print("cluster members:", len(members))

need_by_prefix = defaultdict(set)
for mem in members:
    z = parse_member(mem)
    if z is None:
        print("  UNPARSEABLE member:", mem)
        continue
    arm, raw_gid, genome_prefix, local_pid = z
    need_by_prefix[genome_prefix].add(local_pid)

faa_by_prefix = {}
for path in glob.glob(os.path.join(PROTEOMES, "*.faa")):
    stem = os.path.basename(path)[:-4]
    if stem in need_by_prefix:
        faa_by_prefix[stem] = path

missing = set(need_by_prefix) - set(faa_by_prefix)
if missing:
    sys.exit("STOP: missing proteome files: %s" % sorted(missing))

info = {}
for genome_prefix, wanted_local in need_by_prefix.items():
    path = faa_by_prefix[genome_prefix]
    cur = None
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                local_pid = header.split()[0]
                if local_pid not in wanted_local:
                    cur = None
                    continue
                m = re.search(r"partial=([01]{2})", header)
                partial = m.group(1) if m else "NA"
                cur = genome_prefix + "|" + local_pid
                info[cur] = {"partial": partial, "length": 0}
            elif cur is not None:
                info[cur]["length"] += len(line.strip())

print("")
print("member details (arm, genome, protein, length aa, partial):")
for mem in sorted(members):
    z = parse_member(mem)
    arm, raw_gid, genome_prefix, local_pid = z
    d = info.get(mem, {"partial": "??", "length": -1})
    print("  %-14s %-28s %-30s %5d  %s"
          % (arm, core_genome_id(raw_gid), local_pid,
             d["length"], d["partial"]))

print("")
print("raw DIAMOND hit rows for these members (before qcov filter):")
memset = set(members)
found = 0
with open(HITS) as fh:
    for line in fh:
        q = line.split("\t", 1)[0]
        if q in memset:
            found += 1
            f = line.rstrip("\n").split("\t")
            qcov = 100.0 * float(f[3]) / float(f[4])
            print("  q=%s s=%s pident=%s alen=%s qlen=%s qcov=%.1f e=%s"
                  % (f[0], f[1], f[2], f[3], f[4], qcov, f[6]))
if found == 0:
    print("  none: the amphibian members produced zero alignments at e<=1e-5")

print("")
print("READ IT THIS WAY:")
print("If the rep genome is in the 25 references, the family exists in a")
print("reference by construction and this cluster is not novel content.")
print("Short lengths and partial flags on the amphibian members explain a")
print("missed or low-coverage alignment. Filtered rows shown above mean the")
print("qcov threshold, not absence, removed the hit.")

# ANGELAKISELLA_NOHIT_CLUSTER_CHECK_V1
