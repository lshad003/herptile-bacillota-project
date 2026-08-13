#!/usr/bin/env python3
# Specific and shared clusters are classified by reference homology with shared clusters as the positive control.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/angelakisella_diamond_classify.py
# Output: results/angelakisella_specific_vs_refs_diamond.tsv
"""
Classifies the 67 amphibian-specific and 73 shared Angelakisella MMseqs
clusters by whether their amphibian member proteins hit the 25 GTDB
reference Angelakisella proteomes in the DIAMOND search run as SLURM job
27427721.

Shared clusters are the positive control and must be near 100% hit rate.

Reads:
  results/angelakisella_gene_content.tsv
  work/all_arms_pangenome/mmseqs/all_arms_clusters.tsv
  work/rep_tree/figure_tree_metadata_genus.tsv
  work/angelakisella_diamond/hits.tsv

Writes:
  results/angelakisella_specific_vs_refs_diamond.tsv
Refuses to overwrite.
"""

import csv
import os
import sys
from collections import Counter, defaultdict

BASE = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"

GC = os.path.join(BASE, "results/angelakisella_gene_content.tsv")
CLU = os.path.join(BASE, "work/all_arms_pangenome/mmseqs/all_arms_clusters.tsv")
META = os.path.join(BASE, "work/rep_tree/figure_tree_metadata_genus.tsv")
HITS = os.path.join(BASE, "work/angelakisella_diamond/hits.tsv")
OUT = os.path.join(BASE, "results/angelakisella_specific_vs_refs_diamond.tsv")

GENUS = "Angelakisella"
KEEP_ARMS = {"herptile", "ehi_amphibian", "gtdb_ref"}
AMPH_ARMS = {"herptile", "ehi_amphibian"}

MIN_QCOV = 50.0

if os.path.exists(OUT):
    sys.exit("STOP: output exists, refusing to overwrite:\n" + OUT)

for f in (GC, CLU, META, HITS):
    if not os.path.isfile(f) or os.path.getsize(f) == 0:
        sys.exit("STOP: missing or empty input:\n" + f)


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

n_amph = sum(1 for arm, g in allowed if arm in AMPH_ARMS)
n_ref = sum(1 for arm, g in allowed if arm == "gtdb_ref")
if n_amph != 31 or n_ref != 25:
    sys.exit("STOP: expected 31 amphibian + 25 reference genomes.")


spec = set()
shared = set()
with open(GC) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        fa = float(r["prev_amphibian"])
        fr = float(r["prev_reference"])
        if fa >= 0.90 and fr <= 0.10:
            spec.add(r["cluster"])
        elif fa >= 0.90 and fr >= 0.90:
            shared.add(r["cluster"])

if len(spec) != 67:
    sys.exit("STOP: expected 67 amphibian-specific clusters.")

want = spec | shared
cls = {}
for c in spec:
    cls[c] = "specific"
for c in shared:
    cls[c] = "shared"


def parse_member(pid):
    if "|" not in pid:
        return None
    genome_prefix, local_pid = pid.split("|", 1)
    if "__" not in genome_prefix:
        return None
    arm, raw_gid = genome_prefix.split("__", 1)
    if arm not in KEEP_ARMS:
        return None
    return arm, core_genome_id(raw_gid), genome_prefix, local_pid


query_cluster = {}
n = 0
with open(CLU) as fh:
    for line in fh:
        rep, mem = line.rstrip("\n").split("\t")
        n += 1
        if rep not in want:
            continue
        z = parse_member(mem)
        if z is None:
            continue
        arm, gid, genome_prefix, local_pid = z
        if arm not in AMPH_ARMS or (arm, gid) not in allowed:
            continue
        query_cluster[genome_prefix + "|" + local_pid] = rep

print("cluster rows scanned:", n)
print("amphibian query proteins mapped:", len(query_cluster))
if len(query_cluster) != 4510:
    sys.exit("STOP: expected 4510 query proteins, got %d." % len(query_cluster))

covered = set(query_cluster.values())
if want - covered:
    sys.exit("STOP: %d focal clusters lost their members." % len(want - covered))


best = {}
nhit_rows = 0
with open(HITS) as fh:
    for line in fh:
        q, s, pident, alen, qlen, slen, ev, bits = line.rstrip("\n").split("\t")
        nhit_rows += 1
        c = query_cluster.get(q)
        if c is None:
            sys.exit("STOP: hit query not in focal set:\n" + q)
        qcov = 100.0 * float(alen) / float(qlen)
        if qcov < MIN_QCOV:
            continue
        row = (float(bits), float(pident), qcov, float(ev), q, s)
        if c not in best or row[0] > best[c][0]:
            best[c] = row

print("hit rows read:", nhit_rows)
if nhit_rows != 22004:
    sys.exit("STOP: expected 22004 hit rows from job 27427721.")


summary = Counter()
with open(OUT, "w") as out:
    out.write(
        "cluster\tclass\thit_in_reference\tbest_pident\t"
        "best_qcov\tbest_evalue\tbest_query\tbest_subject\n"
    )
    for c in sorted(want):
        if c in best:
            bits, pident, qcov, ev, q, s = best[c]
            out.write(
                "%s\t%s\tyes\t%.1f\t%.1f\t%.2g\t%s\t%s\n"
                % (c, cls[c], pident, qcov, ev, q, s)
            )
            summary[(cls[c], "hit")] += 1
        else:
            out.write("%s\t%s\tno\t\t\t\t\t\n" % (c, cls[c]))
            summary[(cls[c], "no_hit")] += 1

print("")
print("thresholds: evalue <= 1e-5 (in search), query coverage >= %.0f%%"
      % MIN_QCOV)
for k in ("shared", "specific"):
    hit = summary[(k, "hit")]
    no = summary[(k, "no_hit")]
    tot = hit + no
    print(
        "%s clusters: %d/%d with a reference hit (%.1f%%), %d without"
        % (k, hit, tot, 100.0 * hit / tot, no)
    )

pidents = []
with open(OUT) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["class"] == "specific" and r["hit_in_reference"] == "yes":
            pidents.append(float(r["best_pident"]))
if pidents:
    pidents.sort()
    print(
        "specific clusters with hits: identity min %.1f, median %.1f, max %.1f"
        % (pidents[0], pidents[len(pidents) // 2], pidents[-1])
    )

print("")
print("wrote:", OUT)
print("")
print("READ IT THIS WAY:")
print("Shared clusters must be near 100% hit rate. That is the positive")
print("control. If it is low, the search is broken and nothing else counts.")
print("Specific clusters WITH a hit are divergent versions of families the")
print("references carry; MMseqs split them, they are not novel gene content.")
print("Specific clusters WITHOUT a hit are candidates for genuinely")
print("amphibian-lineage-specific families, pending eggNOG on the survivors.")

# ANGELAKISELLA_DIAMOND_CLASSIFY_V1
