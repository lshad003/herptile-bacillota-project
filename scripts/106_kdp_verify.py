#!/usr/bin/env python3
# The kdpAB locus verified from member annotations and contig adjacency.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/angelakisella_kdp_verify.py
# Output: printed report
"""
Verifies the Kdp interpretation of the amphibian-private OGs COG2060 and
COG2216, and reports the full annotation of all 8 amphibian-private OGs
from actual member proteins rather than OG-level description transfer.

For every amphibian-clade genome (21 herptile + 10 EHI newt):
  1. Finds all proteins assigned each of the 8 private OGs
     (Bacteria-level, taxid 2) in its emapper annotations.
  2. Reports KEGG_ko, Preferred_name, Description per protein.
  3. For COG2060/COG2216 members, tests contig adjacency: same contig
     and gene-index distance <= 3 counts as operon-consistent.
Also checks which of the 26 non-amphibian Angelakisella genomes carry
each private OG (expected: <= 2, by the 0.10 threshold), naming them.

KdpA is K01546, KdpB is K01547 (kdpC K01548, kdpD K07646, kdpE K07667).

Read-only; prints only.
"""

import os
import re
import sys
from collections import defaultdict

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
EGG = B + "/results/eggnog"
EGG_REFS = B + "/results/eggnog_refs"
M_G = B + "/results/angelakisella_matrix_genomes_v2.tsv"

PRIVATE = ["30ET7", "32UD1", "33ABJ", "COG2060", "COG2216",
           "COG3853", "COG3882", "COG4915"]
KDP = {"COG2060", "COG2216"}
KDP_KOS = {"ko:K01546": "kdpA", "ko:K01547": "kdpB", "ko:K01548": "kdpC",
           "ko:K07646": "kdpD", "ko:K07667": "kdpE"}

genomes = []
with open(M_G) as fh:
    fh.readline()
    for line in fh:
        g, cat, grp = line.rstrip("\n").split("\t")
        genomes.append((g, cat, grp))

amph = [(g, c) for g, c, gr in genomes if gr == "amphibian"]
nonamph = [(g, c) for g, c, gr in genomes if gr == "non_amphibian"]
print("amphibian %d, non_amphibian %d" % (len(amph), len(nonamph)))
if len(amph) != 31 or len(nonamph) != 26:
    sys.exit("STOP: expected 31 + 26.")


def ann_file(key, cat):
    arm, gid = key.split("__", 1)
    d = EGG_REFS if arm == "gtdb_ref" else EGG
    return os.path.join(d, gid + ".emapper.annotations")


def scan(path, want_ogs):
    """protein -> (set of matched private OGs at taxid 2, ko, name, desc)"""
    hits = []
    hdr = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#query"):
                hdr = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if hdr is None or len(f) < len(hdr):
                continue
            row = dict(zip(hdr, f))
            matched = set()
            for og in row.get("eggNOG_OGs", "").split(","):
                og = og.strip()
                if "@" not in og:
                    continue
                unit, tax = og.split("@", 1)
                if tax.split("|", 1)[0] == "2" and unit in want_ogs:
                    matched.add(unit)
            if matched:
                hits.append((row["query"], matched,
                             row.get("KEGG_ko", "-"),
                             row.get("Preferred_name", "-"),
                             row.get("Description", "-")[:80]))
    return hits


def gene_pos(pid):
    """k141_274340_2 -> (contig, index). Rsplit is safe here because the
    trailing _N is Prodigal's gene index; contig names keep their own
    underscores intact on the left."""
    m = re.match(r"^(.*)_(\d+)$", pid)
    if m is None:
        return None, None
    return m.group(1), int(m.group(2))


ko_counts = defaultdict(lambda: defaultdict(int))
adjacency = []
per_og_examples = defaultdict(list)

for key, cat in amph:
    path = ann_file(key, cat)
    if not os.path.isfile(path):
        sys.exit("STOP: missing annotations: " + path)
    hits = scan(path, set(PRIVATE))
    kdp_pos = []
    for pid, ogs, ko, name, desc in hits:
        for og in ogs:
            ko_counts[og][ko] += 1
            if len(per_og_examples[og]) < 3:
                per_og_examples[og].append(
                    "%s %s ko=%s name=%s | %s" % (key, pid, ko, name, desc))
        if ogs & KDP:
            contig, idx = gene_pos(pid)
            kdp_pos.append((contig, idx, sorted(ogs & KDP), ko))
    both = {}
    for contig, idx, ogs, ko in kdp_pos:
        both.setdefault(contig, []).append((idx, ogs, ko))
    for contig, items in both.items():
        units = set()
        for idx, ogs, ko in items:
            units.update(ogs)
        if units >= KDP:
            idxs = sorted(i for i, o, k in items)
            span = max(idxs) - min(idxs)
            adjacency.append((key, contig, idxs, span))

print("")
print("=" * 70)
print("KO composition of each amphibian-private OG (across 31 genomes):")
for og in PRIVATE:
    print("")
    print(og + ":")
    for ko, c in sorted(ko_counts[og].items(), key=lambda z: -z[1]):
        label = ""
        for k, g in KDP_KOS.items():
            if k in ko:
                label = "  <-- %s" % g
        print("  %3d  %s%s" % (c, ko, label))
    for ex in per_og_examples[og]:
        print("   e.g. " + ex)

print("")
print("=" * 70)
print("COG2060 + COG2216 contig adjacency (operon test):")
if adjacency:
    for key, contig, idxs, span in adjacency:
        verdict = "ADJACENT" if span <= 3 else "same contig, span %d" % span
        print("  %-40s %-22s genes %s  %s"
              % (key, contig, idxs, verdict))
    n_adj = sum(1 for _, _, _, s in adjacency if s <= 3)
    print("genomes with both OGs on one contig: %d; adjacent (<=3): %d"
          % (len(adjacency), n_adj))
else:
    print("  no genome has both OGs on a single contig")
print("NOTE: absence of same-contig co-occurrence in MAGs can reflect")
print("assembly fragmentation; presence is the informative direction.")

print("")
print("=" * 70)
print("Non-amphibian carriers of each private OG (threshold allows <= 2):")
for og in PRIVATE:
    carriers = []
    for key, cat in nonamph:
        path = ann_file(key, cat)
        if not os.path.isfile(path):
            sys.exit("STOP: missing annotations: " + path)
        if scan(path, {og}):
            carriers.append(key)
    tag = ""
    if og == "COG3882" and not any("youngblut" in c for c in carriers):
        tag = "  (Youngblut genome also lacks FkbH: confirmed)"
    print("  %-10s %d carriers%s" % (og, len(carriers), tag))
    for c in carriers:
        print("      " + c)

# ANGELAKISELLA_KDP_VERIFY_V1
