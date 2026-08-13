#!/usr/bin/env python3
# Private units per clade counted in both currencies.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/angelakisella_block_privacy_stats.py
# Output: results/angelakisella_block_privacy.tsv
"""
Quantifies the panel contrast in the neighborhood figure: for the
amphibian clade and every reference subclade of >= MIN_BLOCK genomes,
counts private units (present in >= 0.90 of the block, <= 0.10 of the
non-context remainder) in BOTH currencies, and reports the MMseqs/OG
privacy ratio per block.

Reading it: if the amphibian block's ratio is similar to the reference
subclades' ratios, apparent clade-specific content collapses uniformly
under orthology correction and gene content tracks divergence in every
direction. If the amphibian ratio is much lower (more OG-private
content survives), the amphibian clade retains disproportionate
functional distinctiveness, which is the reportable signature.

Reference subclades are defined from the pruned tree: maximal clades
containing only non_amphibian genomes, of size >= MIN_BLOCK.
Context (Heteroruminococcus) genomes are excluded from both blocks and
background.

Writes results/angelakisella_block_privacy.tsv (refuses overwrite).
Also lists the OG-private units of the amphibian block with their IDs
so the FkbH family can be located among them.
"""

import os
import sys
import numpy as np
import dendropy

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
NWK = B + "/work/rep_tree/figure_tree.nwk"
M_OG = B + "/results/angelakisella_matrix_og_v2.tsv"
M_MM = B + "/results/angelakisella_matrix_mmseqs_v2.tsv"
M_G = B + "/results/angelakisella_matrix_genomes_v2.tsv"
OUT = B + "/results/angelakisella_block_privacy.tsv"

MIN_BLOCK = 4
P_IN = 0.90
P_OUT = 0.10

if os.path.exists(OUT):
    sys.exit("STOP: output exists, refusing to overwrite:\n" + OUT)
for p in (NWK, M_OG, M_MM, M_G):
    if not os.path.exists(p):
        sys.exit("STOP: missing input:\n" + p)


def read_matrix(path):
    with open(path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")[1:]
        units, rows = [], []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            units.append(f[0])
            rows.append([int(x) for x in f[1:]])
    return cols, units, np.array(rows, dtype=np.int8)


cols_og, units_og, M1 = read_matrix(M_OG)
cols_mm, units_mm, M2 = read_matrix(M_MM)
if cols_og != cols_mm:
    sys.exit("STOP: matrix columns disagree.")
cols = cols_og

grp = {}
with open(M_G) as fh:
    fh.readline()
    for line in fh:
        g, cat, gr = line.rstrip("\n").split("\t")
        grp[g] = gr

col_i = {g: i for i, g in enumerate(cols)}
amph = [g for g in cols if grp[g] == "amphibian"]
nonamph = [g for g in cols if grp[g] == "non_amphibian"]
context = [g for g in cols if grp[g] == "context"]
print("amphibian %d, non_amphibian %d, context %d"
      % (len(amph), len(nonamph), len(context)))
if (len(amph), len(nonamph), len(context)) != (31, 26, 4):
    sys.exit("STOP: expected 31/26/4.")


def tip_key(label):
    f = label.split("|")
    if len(f) < 2:
        return None
    gid = f[-1]
    if gid.startswith(("RS_", "GB_")):
        gid = gid[3:]
    return f[0] + "__" + gid


t = dendropy.Tree.get(path=NWK, schema="newick", preserve_underscores=True)
keep = set(cols)
allb = [x.label for x in t.taxon_namespace if x.label is not None]
t.prune_taxa_with_labels([b for b in allb if tip_key(b) not in keep])
if hasattr(t, "purge_taxon_namespace"):
    t.purge_taxon_namespace()
t.suppress_unifurcations()
if len(t.leaf_nodes()) != 61:
    sys.exit("STOP: pruned tree tip count != 61.")

# maximal all-non_amphibian clades of size >= MIN_BLOCK
blocks = []
visited = set()
for nd in t.preorder_node_iter():
    if id(nd) in visited:
        continue
    tips = [tip_key(lf.taxon.label) for lf in nd.leaf_iter()]
    if all(grp.get(k) == "non_amphibian" for k in tips):
        if len(tips) >= MIN_BLOCK:
            blocks.append(("ref_subclade_%d" % (len(blocks) + 1), tips))
        for x in nd.preorder_iter():
            visited.add(id(x))

blocks.insert(0, ("amphibian_clade", amph))
print("blocks: %s" % ", ".join("%s(%d)" % (n, len(g)) for n, g in blocks))

background_pool = set(amph) | set(nonamph)


def privacy(M, units, block_genomes):
    bi = [col_i[g] for g in block_genomes]
    oi = [col_i[g] for g in background_pool - set(block_genomes)]
    inprev = M[:, bi].mean(axis=1)
    outprev = M[:, oi].mean(axis=1)
    mask = (inprev >= P_IN) & (outprev <= P_OUT)
    return [units[i] for i in np.where(mask)[0]]


rows_out = []
amph_og_private = None
for name, gset in blocks:
    p_og = privacy(M1, units_og, gset)
    p_mm = privacy(M2, units_mm, gset)
    ratio = (len(p_mm) / len(p_og)) if p_og else float("inf")
    rows_out.append((name, len(gset), len(p_og), len(p_mm), ratio))
    if name == "amphibian_clade":
        amph_og_private = p_og

print("")
print("%-20s %6s %10s %12s %10s"
      % ("block", "n", "OG_priv", "MMseqs_priv", "MM/OG"))
with open(OUT, "w") as out:
    out.write("block\tn_genomes\tog_private\tmmseqs_private\tratio\n")
    for name, n, og, mm, r in rows_out:
        rtxt = "inf" if r == float("inf") else "%.2f" % r
        print("%-20s %6d %10d %12d %10s" % (name, n, og, mm, rtxt))
        out.write("%s\t%d\t%d\t%d\t%s\n" % (name, n, og, mm, rtxt))

print("")
print("wrote:", OUT)

if amph_og_private is not None:
    print("")
    print("amphibian-clade OG-private units (%d):" % len(amph_og_private))
    for u in amph_og_private:
        print("  ", u)
    print("")
    print("To find the FkbH family among these: its herptile members were")
    print("annotated COG3882; check whether COG3882 appears above.")

print("")
print("Thresholds: private = prevalence >= %.2f in block, <= %.2f in the"
      % (P_IN, P_OUT))
print("other Angelakisella genomes. Context genomes excluded throughout.")

# ANGELAKISELLA_BLOCK_PRIVACY_V1
