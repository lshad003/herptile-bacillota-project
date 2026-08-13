#!/usr/bin/env python3
# Amphibian member proteins of specific and shared clusters are staged as DIAMOND queries with the 25 reference proteomes as the database.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/angelakisella_specific_vs_refs_diamond.py
# Output: work/angelakisella_diamond/amphibian_queries.faa, work/angelakisella_diamond/reference_proteomes.faa
"""
Tests whether the 67 near-fixed amphibian-specific Angelakisella MMseqs
clusters are genuinely absent from the 25 GTDB reference Angelakisella
genomes, or are divergent versions of gene families the references carry.

Method:
  1. Rebuild the focal genome set and the specific/shared cluster sets
     exactly as in angelakisella_cluster_quality.py.
  2. Extract amphibian member proteins of BOTH cluster classes as queries.
     Shared clusters are the positive control: their amphibian members
     must hit the reference proteomes.
  3. Build a DIAMOND database from the 25 reference proteomes only.
  4. blastp --very-sensitive, classify each cluster by its best hit.

Writes:
  results/angelakisella_specific_vs_refs_diamond.tsv  (per-cluster)
Refuses to overwrite existing outputs.
"""

import csv
import glob
import os
import shutil
import subprocess
import sys
from collections import Counter, defaultdict

BASE = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"

GC = os.path.join(BASE, "results/angelakisella_gene_content.tsv")
CLU = os.path.join(BASE, "work/all_arms_pangenome/mmseqs/all_arms_clusters.tsv")
PROTEOMES = os.path.join(BASE, "work/all_arms_pangenome/proteomes")
META = os.path.join(BASE, "work/rep_tree/figure_tree_metadata_genus.tsv")

WORK = os.path.join(BASE, "work/angelakisella_diamond")
OUT = os.path.join(BASE, "results/angelakisella_specific_vs_refs_diamond.tsv")

GENUS = "Angelakisella"
KEEP_ARMS = {"herptile", "ehi_amphibian", "gtdb_ref"}
AMPH_ARMS = {"herptile", "ehi_amphibian"}

EVALUE = 1e-5
MIN_QCOV = 50.0

if os.path.exists(OUT):
    sys.exit("STOP: output exists, refusing to overwrite:\n" + OUT)

diamond = shutil.which("diamond")
if diamond is None:
    sys.exit(
        "STOP: diamond not found in PATH.\n"
        "Run: module load diamond\n"
        "then rerun this script."
    )
print("diamond binary:", diamond)

os.makedirs(WORK, exist_ok=True)


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


# ------------------------------------------------------------
# 1. Focal genomes
# ------------------------------------------------------------

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
print("genomes: %d amphibian, %d reference" % (n_amph, n_ref))
if n_amph != 31 or n_ref != 25:
    sys.exit("STOP: expected 31 amphibian + 25 reference genomes.")


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


# ------------------------------------------------------------
# 2. Cluster classes
# ------------------------------------------------------------

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

print("specific clusters: %d, shared clusters: %d" % (len(spec), len(shared)))
if len(spec) != 67:
    sys.exit("STOP: expected 67 amphibian-specific clusters.")

want = spec | shared
cls = {}
for c in spec:
    cls[c] = "specific"
for c in shared:
    cls[c] = "shared"


# ------------------------------------------------------------
# 3. Amphibian member proteins of those clusters
# ------------------------------------------------------------

query_pids = {}
need_by_prefix = defaultdict(set)
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
        query_pids[genome_prefix + "|" + local_pid] = rep
        need_by_prefix[genome_prefix].add(local_pid)
        if n % 2000000 == 0:
            print("  scanned %d cluster rows" % n)

print("cluster rows scanned:", n)
print("amphibian query proteins:", len(query_pids))
covered = set(query_pids.values())
missing = want - covered
if missing:
    sys.exit("STOP: %d focal clusters have no amphibian members." % len(missing))


# ------------------------------------------------------------
# 4. Write query FASTA from original proteomes
# ------------------------------------------------------------

faa_by_prefix = {}
ref_faas = []
for path in glob.glob(os.path.join(PROTEOMES, "*.faa")):
    stem = os.path.basename(path)[:-4]
    if stem in need_by_prefix:
        faa_by_prefix[stem] = path
    z = None
    if "__" in stem:
        arm, raw_gid = stem.split("__", 1)
        if arm == "gtdb_ref" and (arm, core_genome_id(raw_gid)) in allowed:
            ref_faas.append(path)

print("query proteome files found: %d of %d"
      % (len(faa_by_prefix), len(need_by_prefix)))
print("reference proteome files found: %d of 25" % len(ref_faas))
if len(faa_by_prefix) != len(need_by_prefix):
    sys.exit("STOP: missing query proteome files.")
if len(ref_faas) != 25:
    sys.exit("STOP: expected 25 reference proteome files.")

query_faa = os.path.join(WORK, "amphibian_queries.faa")
nq = 0
with open(query_faa, "w") as out:
    for genome_prefix in sorted(need_by_prefix):
        wanted_local = need_by_prefix[genome_prefix]
        keep = False
        with open(faa_by_prefix[genome_prefix], errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    local_pid = line[1:].split()[0]
                    keep = local_pid in wanted_local
                    if keep:
                        out.write(">" + genome_prefix + "|" + local_pid + "\n")
                        nq += 1
                elif keep:
                    out.write(line)
print("query proteins written:", nq)
if nq != len(query_pids):
    sys.exit("STOP: query FASTA count mismatch.")

ref_faa = os.path.join(WORK, "reference_proteomes.faa")
nr = 0
with open(ref_faa, "w") as out:
    for path in sorted(ref_faas):
        stem = os.path.basename(path)[:-4]
        with open(path, errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    local_pid = line[1:].split()[0]
                    out.write(">" + stem + "|" + local_pid + "\n")
                    nr += 1
                else:
                    out.write(line)
print("reference proteins written:", nr)


# ------------------------------------------------------------
# 5. DIAMOND
# ------------------------------------------------------------

db = os.path.join(WORK, "ref_db")
hits = os.path.join(WORK, "hits.tsv")

subprocess.run(
    [diamond, "makedb", "--in", ref_faa, "-d", db, "--quiet"],
    check=True,
)
subprocess.run(
    [
        diamond, "blastp", "--very-sensitive", "--quiet",
        "-q", query_faa, "-d", db, "-o", hits,
        "-e", str(EVALUE), "-k", "5",
        "--outfmt", "6",
        "qseqid", "sseqid", "pident", "length",
        "qlen", "slen", "evalue", "bitscore",
    ],
    check=True,
)
print("diamond finished")


# ------------------------------------------------------------
# 6. Best hit per cluster
# ------------------------------------------------------------

best = {}
with open(hits) as fh:
    for line in fh:
        q, s, pident, alen, qlen, slen, ev, bits = line.rstrip("\n").split("\t")
        c = query_pids.get(q)
        if c is None:
            continue
        qcov = 100.0 * float(alen) / float(qlen)
        if qcov < MIN_QCOV:
            continue
        row = (float(bits), float(pident), qcov, float(ev), q, s)
        if c not in best or row[0] > best[c][0]:
            best[c] = row

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
print("thresholds: evalue <= %g, query coverage >= %.0f%%" % (EVALUE, MIN_QCOV))
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
print("READ IT THIS WAY:")
print("Shared clusters must be near 100%% hit rate. That is the positive")
print("control. If it is low, the search is broken and nothing else counts.")
print("Specific clusters WITH a hit are divergent versions of families the")
print("references carry; MMseqs split them, they are not novel gene content.")
print("Specific clusters WITHOUT a hit are candidates for genuinely")
print("amphibian-lineage-specific families, pending eggNOG on the survivors.")

# ANGELAKISELLA_SPECIFIC_VS_REFS_DIAMOND_V1
