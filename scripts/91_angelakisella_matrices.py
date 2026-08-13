#!/usr/bin/env python3
# The 61-genome Angelakisella neighborhood is defined from the tree and OG and MMseqs presence matrices are built.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/angelakisella_build_matrices_v2.py
# Output: results/angelakisella_matrix_og_v2.tsv, results/angelakisella_matrix_mmseqs_v2.tsv, results/angelakisella_matrix_genomes_v2.tsv
"""
V2: builds presence matrices for the 61-tip Angelakisella neighborhood
figure (supersedes the unbuilt 56-genome V1 plan).

Genome set = every tip under the neighborhood node (MRCA of all
Angelakisella tips in the figure tree, plus 1 parent level):
  21 wild herptile + 10 EHI newt   (amphibian clade)
  25 GTDB reference + 1 Youngblut  (non-amphibian Angelakisella)
  4 Heteroruminococcus references  (context)

Currencies:
  1. eggNOG Bacteria-level OGs (taxid 2), per-genome annotations
  2. raw MMseqs clusters, restricted to the 61 genomes (divergence
     comparison panel only)

Units in fewer than MIN_GENOMES of the 61 are dropped.

GATE: refuses to run unless all 61 annotation files exist and are
nonempty (waits for emapper arrays 27428475 and 27430034).

NOTE: the Youngblut genome ID contains "__"; genome IDs are treated as
literal strings, arm prefixes are split on the FIRST "__" only.

Writes:
  results/angelakisella_matrix_og_v2.tsv
  results/angelakisella_matrix_mmseqs_v2.tsv
  results/angelakisella_matrix_genomes_v2.tsv  (genome, catalog, group)
Refuses to overwrite existing outputs.
"""

import csv
import os
import sys
from collections import defaultdict

import dendropy

BASE = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
NWK = os.path.join(BASE, "work/rep_tree/figure_tree.nwk")
META = os.path.join(BASE, "work/rep_tree/figure_tree_metadata_genus.tsv")
EGG = os.path.join(BASE, "results/eggnog")
EGG_REFS = os.path.join(BASE, "results/eggnog_refs")
CLU = os.path.join(BASE, "work/all_arms_pangenome/mmseqs/all_arms_clusters.tsv")

OUT_OG = os.path.join(BASE, "results/angelakisella_matrix_og_v2.tsv")
OUT_MM = os.path.join(BASE, "results/angelakisella_matrix_mmseqs_v2.tsv")
OUT_G = os.path.join(BASE, "results/angelakisella_matrix_genomes_v2.tsv")

GENUS = "Angelakisella"
LEVELS_UP = 1
MIN_GENOMES = 3
AMPH_ARMS = {"herptile", "ehi_amphibian"}
ALL_ARMS = {"herptile", "ehi_amphibian", "gtdb_ref", "youngblut", "ehi"}

for p in (OUT_OG, OUT_MM, OUT_G):
    if os.path.exists(p):
        sys.exit("STOP: output exists, refusing to overwrite:\n" + p)


def norm_gid(raw, arm):
    x = raw.strip()
    if x.startswith(arm + "__"):
        x = x[len(arm) + 2:]
    if x.startswith(arm + "|"):
        x = x[len(arm) + 1:]
    if x.startswith(arm + "|"):
        x = x[len(arm) + 1:]
    if x.startswith(("RS_", "GB_")):
        x = x[3:]
    return x


genus_of = {}
arm_of = {}
with open(META) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        arm = r["arm"].strip()
        k = arm + "__" + norm_gid(r["genome"], arm)
        genus_of[k] = r["genus"].strip()
        arm_of[k] = arm

focal = {k for k, g in genus_of.items() if g == GENUS}
print("%s genomes in metadata: %d" % (GENUS, len(focal)))
if len(focal) != 57:
    sys.exit("STOP: expected 57 Angelakisella genomes in tree metadata.")


def tip_key(label):
    f = label.split("|")
    if len(f) < 2:
        return None
    gid = f[-1]
    if gid.startswith(("RS_", "GB_")):
        gid = gid[3:]
    return f[0] + "__" + gid


t = dendropy.Tree.get(path=NWK, schema="newick", preserve_underscores=True)
fl = [lf for lf in t.leaf_nodes()
      if lf.taxon and tip_key(lf.taxon.label) in focal]
print("focal tips in tree: %d" % len(fl))
if len(fl) != 57:
    sys.exit("STOP: expected 57 focal tips in the tree.")

fid = set(id(x) for x in fl)
node = fl[0]
while node is not None:
    under = set(id(x) for x in node.leaf_iter())
    if fid <= under:
        break
    node = node.parent_node
for _ in range(LEVELS_UP):
    if node.parent_node is not None:
        node = node.parent_node

genomes = []
for lf in node.leaf_iter():
    k = tip_key(lf.taxon.label)
    genomes.append(k)
print("neighborhood genomes: %d" % len(genomes))
if len(genomes) != 61:
    sys.exit("STOP: expected 61 neighborhood tips, got %d." % len(genomes))


def group_of(k):
    arm = arm_of.get(k, "?")
    g = genus_of.get(k, "?")
    if g == GENUS and arm in AMPH_ARMS:
        return "amphibian"
    if g == GENUS:
        return "non_amphibian"
    return "context"


n_a = sum(1 for k in genomes if group_of(k) == "amphibian")
n_n = sum(1 for k in genomes if group_of(k) == "non_amphibian")
n_c = sum(1 for k in genomes if group_of(k) == "context")
print("groups: %d amphibian, %d non_amphibian, %d context" % (n_a, n_n, n_c))
if (n_a, n_n, n_c) != (31, 26, 4):
    sys.exit("STOP: expected 31 + 26 + 4.")

# ---------------------------------------------------------- gate
ann_path = {}
missing = []
for k in genomes:
    arm, gid = k.split("__", 1)
    d = EGG_REFS if arm == "gtdb_ref" else EGG
    p = os.path.join(d, gid + ".emapper.annotations")
    if os.path.isfile(p) and os.path.getsize(p) > 0:
        ann_path[k] = p
    else:
        missing.append(k)

if missing:
    print("")
    for k in missing:
        print("MISSING:", k)
    sys.exit("STOP: %d/61 annotation files missing; wait for arrays "
             "27428475 and 27430034 and rerun." % len(missing))

print("gate passed: 61/61 annotation files present")

# ---------------------------------------------------------- OG matrix
og_presence = defaultdict(set)
per_genome = {}
for k in sorted(genomes):
    n_rows = 0
    hdr = None
    with open(ann_path[k]) as fh:
        for line in fh:
            if line.startswith("#query"):
                hdr = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if hdr is None or len(f) < len(hdr):
                continue
            n_rows += 1
            row = dict(zip(hdr, f))
            for og in row.get("eggNOG_OGs", "").split(","):
                og = og.strip()
                if "@" not in og:
                    continue
                unit, tax = og.split("@", 1)
                if tax.split("|", 1)[0] == "2":
                    og_presence[unit].add(k)
                    break
    per_genome[k] = n_rows

print("OG units before prevalence filter: %d" % len(og_presence))
zero = [k for k, v in per_genome.items() if v == 0]
if zero:
    for k in zero:
        print("EMPTY ANNOTATIONS:", k)
    sys.exit("STOP: genomes with zero annotation rows.")

# ---------------------------------------------------------- MMseqs matrix
gset = set(genomes)


def parse_member(pid):
    if "|" not in pid:
        return None
    genome_prefix, local_pid = pid.split("|", 1)
    if "__" not in genome_prefix:
        return None
    arm, raw_gid = genome_prefix.split("__", 1)
    if arm not in ALL_ARMS:
        return None
    if raw_gid.startswith(("RS_", "GB_")):
        raw_gid = raw_gid[3:]
    return arm + "__" + raw_gid


mm_presence = defaultdict(set)
n = 0
with open(CLU) as fh:
    for line in fh:
        rep, mem = line.rstrip("\n").split("\t")
        n += 1
        k = parse_member(mem)
        if k is None or k not in gset:
            continue
        mm_presence[rep].add(k)
        if n % 2000000 == 0:
            print("  scanned %d cluster rows" % n)

print("cluster rows scanned: %d" % n)
print("MMseqs units before prevalence filter: %d" % len(mm_presence))

genomes_missing_mm = gset - set().union(*mm_presence.values()) \
    if mm_presence else gset
if genomes_missing_mm:
    for k in sorted(genomes_missing_mm):
        print("NO MMSEQS MEMBERSHIP:", k)
    sys.exit("STOP: genomes absent from the cluster tsv.")


def write_matrix(path, presence, label):
    units = sorted(u for u, gs in presence.items()
                   if len(gs) >= MIN_GENOMES)
    cols = sorted(genomes)
    with open(path, "w") as out:
        out.write("unit\t" + "\t".join(cols) + "\n")
        for u in units:
            gs = presence[u]
            out.write(u + "\t" + "\t".join(
                "1" if c in gs else "0" for c in cols) + "\n")
    print("%s: %d units (>= %d genomes) x %d genomes -> %s"
          % (label, len(units), MIN_GENOMES, len(cols), path))
    return len(units)


n_og = write_matrix(OUT_OG, og_presence, "OG matrix")
n_mm = write_matrix(OUT_MM, mm_presence, "MMseqs matrix")

with open(OUT_G, "w") as out:
    out.write("genome\tcatalog\tgroup\n")
    for k in sorted(genomes):
        out.write("%s\t%s\t%s\n" % (k, arm_of.get(k, "?"), group_of(k)))
print("genome metadata -> %s" % OUT_G)

print("")
print("Unit-count ratio MMseqs/OG: %.2f" % (n_mm / max(n_og, 1)))

# ANGELAKISELLA_BUILD_MATRICES_V2
