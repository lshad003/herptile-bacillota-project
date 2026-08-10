#!/usr/bin/env python3
# Dereplication clusters are audited for arm composition, giving the cross-arm positive control.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/pooled_drep_multiarm_audit.py
# Output: results/pooled_drep_multiarm_audit.tsv
import os, sys, csv
from collections import defaultdict, Counter

BASE = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
CDB  = os.path.join(BASE, "work/pooled_drep_rum/drep_out/data_tables/Cdb.csv")
ARMS = os.path.join(BASE, "work/pooled_drep_rum/genome_arms.tsv")
OUT  = os.path.join(BASE, "results/pooled_drep_multiarm_audit.tsv")

def die(m):
    sys.stderr.write("FATAL: " + m + "\n"); sys.exit(1)

if os.path.exists(OUT):
    die("output exists, refusing to overwrite: " + OUT)
for p in (CDB, ARMS):
    if not os.path.isfile(p):
        die("missing " + p)

arm_of = {}; gid_of = {}; nb = 0
with open(ARMS) as fh:
    head = fh.readline().rstrip("\n").split("\t")
    for c in ("staged_name", "arm", "genome_id"):
        if c not in head:
            die("genome_arms.tsv lacks column " + c)
    si, ai, gi = head.index("staged_name"), head.index("arm"), head.index("genome_id")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if len(f) <= max(si, ai, gi):
            continue
        s = f[si]
        if s in arm_of:
            die("duplicate staged_name in bridge: " + s)
        arm_of[s] = f[ai]; gid_of[s] = f[gi]; nb += 1
if nb != 1957:
    die("expected 1957 bridge rows, found %d" % nb)

EXT = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz",
       ".contigs.fa", ".contigs.fa.gz")

def resolve(name):
    if name in arm_of:
        return name, False
    for e in sorted(EXT, key=len, reverse=True):
        if name.endswith(e) and name[:-len(e)] in arm_of:
            return name[:-len(e)], True
    return None, False

rows = 0; n_strip = 0
clusters = defaultdict(list)
with open(CDB, newline="") as fh:
    rd = csv.DictReader(fh)
    fn = rd.fieldnames or []
    for c in ("genome", "secondary_cluster"):
        if c not in fn:
            die("Cdb.csv lacks '%s'; columns are: %s" % (c, ",".join(fn)))
    for r in rd:
        rows += 1
        key, stripped = resolve(r["genome"])
        if key is None:
            die("Cdb genome absent from bridge: " + r["genome"])
        if stripped:
            n_strip += 1
        clusters[r["secondary_cluster"]].append((key, arm_of[key], gid_of[key]))
if rows != 1895:
    die("expected 1895 Cdb rows, found %d" % rows)

print("bridge rows %d | Cdb rows %d | filename extension stripped on %d" % (nb, rows, n_strip))
print("bridge minus Cdb = %d (expected 62 dropped by -comp 70)" % (nb - rows))
print("")

bridge_arm = Counter(arm_of.values())
cdb_arm = Counter(a for v in clusters.values() for (_, a, _) in v)
print("=== ARM LABELS AS FOUND, genomes not units ===")
for a in sorted(bridge_arm):
    print("  %-18s bridge %5d | clustered %5d" % (a, bridge_arm[a], cdb_arm.get(a, 0)))
print("")

n_clu = len(clusters)
multi = {c: v for c, v in clusters.items() if len(set(a for _, a, _ in v)) > 1}
print("=== CLUSTER COUNTS ===")
print("  clusters %d   (CHATINDEX records 1,515)  %s"
      % (n_clu, "MATCH" if n_clu == 1515 else "MISMATCH, DO NOT USE BELOW"))
print("  single-arm %d" % (n_clu - len(multi)))
print("  multi-arm  %d   (CHATINDEX records 57)   %s"
      % (len(multi), "MATCH" if len(multi) == 57 else "MISMATCH, DO NOT USE BELOW"))
print("")

combos = Counter(tuple(sorted(set(a for _, a, _ in v))) for v in multi.values())
print("=== ARM COMBINATIONS AMONG MULTI-ARM CLUSTERS ===")
for k, n in sorted(combos.items(), key=lambda x: (-x[1], x[0])):
    print("  %3d  %s" % (n, " + ".join(k)))
print("")

print("=== PER ARM: CLUSTERS SHARED WITH A DIFFERENT ARM ===")
clean = []
for a in sorted(bridge_arm):
    tot = sum(1 for v in clusters.values() if any(x[1] == a for x in v))
    sh  = sum(1 for v in multi.values() if any(x[1] == a for x in v))
    print("  %-18s in %5d clusters | shared with another arm %3d" % (a, tot, sh))
    if tot > 0 and sh == 0:
        clean.append(a)
print("")
print("ARMS SHARING NO CLUSTER WITH ANY OTHER ARM:")
if clean:
    for a in clean:
        print("  " + a)
else:
    print("  none")
print("")
print("READ THIS AGAINST THE 3.2 SENTENCE. The claim holds only for arms")
print("listed immediately above. Name the arm as it appears in the table.")
print("")

with open(OUT, "w") as out:
    out.write("secondary_cluster\tn_genomes\tn_arms\tarms\tgenome_ids\n")
    for c in sorted(multi):
        v = multi[c]
        out.write("%s\t%d\t%d\t%s\t%s\n" % (
            c, len(v), len(set(x[1] for x in v)),
            ";".join(sorted(set(x[1] for x in v))),
            ";".join(x[2] for x in v)))
print("WROTE: " + OUT)
# POOLED_DREP_MULTIARM_V1_20260808
