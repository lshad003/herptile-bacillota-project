#!/usr/bin/env python3
# Neighborhood subtree with orthologous-group and cluster presence heatmaps.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/fig_angelakisella_neighborhood_heatmap.py
# Output: results/figures/Figure_angelakisella_neighborhood_heatmap.pdf/.png
"""
Merged Angelakisella neighborhood figure: the 61-tip bac120 subtree
(57 Angelakisella + 4 Heteroruminococcus context) with two aligned
presence heatmaps sharing the tree row order:

  panel 1: eggNOG Bacteria-level OGs (corrected currency)
  panel 2: raw MMseqs clusters (divergence comparison)

Supersedes the separate tree-only neighborhood figure and the 56-genome
tree+heatmap plan.

Reads:
  work/rep_tree/figure_tree.nwk
  results/angelakisella_matrix_og_v2.tsv
  results/angelakisella_matrix_mmseqs_v2.tsv
  results/angelakisella_matrix_genomes_v2.tsv
Writes:
  results/figures/Figure_angelakisella_neighborhood_heatmap.pdf/.png
Refuses to overwrite. Login node fine.
"""

import os
import sys
import numpy as np
import dendropy
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, to_rgba
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
NWK = B + "/work/rep_tree/figure_tree.nwk"
M_OG = B + "/results/angelakisella_matrix_og_v2.tsv"
M_MM = B + "/results/angelakisella_matrix_mmseqs_v2.tsv"
M_G = B + "/results/angelakisella_matrix_genomes_v2.tsv"
OUTD = B + "/results/figures"
STEM = OUTD + "/Figure_angelakisella_neighborhood_heatmap"

CATC = {
    "herptile": to_rgba("#C44E52"),
    "ehi_amphibian": to_rgba("#DD8452"),
    "youngblut": to_rgba("#8172B3"),
    "gtdb_ref": to_rgba("#4C72B0"),
}
CATL = {
    "herptile": "wild herptile",
    "ehi_amphibian": "EHI newt",
    "youngblut": "Youngblut",
    "gtdb_ref": "GTDB reference",
}
GRPC = {
    "amphibian": to_rgba("#C44E52"),
    "non_amphibian": to_rgba("#4C72B0"),
    "context": to_rgba("#BEBEBE"),
}
GRPL = {
    "amphibian": "amphibian Angelakisella",
    "non_amphibian": "other Angelakisella",
    "context": "Heteroruminococcus",
}
PRESENT = to_rgba("#2B4C7E")
ABSENT = to_rgba("#E8E8E8")

for p in (NWK, M_OG, M_MM, M_G):
    if not os.path.exists(p):
        sys.exit("STOP: missing input (run the V2 matrix builder first):\n" + p)
for ext in (".pdf", ".png"):
    if os.path.exists(STEM + ext):
        sys.exit("STOP: output exists, refusing to overwrite:\n" + STEM + ext)
os.makedirs(OUTD, exist_ok=True)


def read_matrix(path):
    with open(path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")[1:]
        units = []
        rows = []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            units.append(f[0])
            rows.append([int(x) for x in f[1:]])
    return cols, units, np.array(rows, dtype=np.int8)


cols_og, units_og, M1 = read_matrix(M_OG)
cols_mm, units_mm, M2 = read_matrix(M_MM)
if cols_og != cols_mm:
    sys.exit("STOP: matrix genome columns disagree.")
cols = cols_og
print("OG matrix: %d units x %d genomes" % (len(units_og), len(cols)))
print("MMseqs matrix: %d units x %d genomes" % (len(units_mm), len(cols)))

meta = {}
with open(M_G) as fh:
    fh.readline()
    for line in fh:
        g, cat, grp = line.rstrip("\n").split("\t")
        meta[g] = (cat, grp)
if set(meta) != set(cols):
    sys.exit("STOP: genome metadata does not match matrix columns.")

# ---------------------------------------------------------------- tree
t = dendropy.Tree.get(path=NWK, schema="newick", preserve_underscores=True)


def tip_key(label):
    f = label.split("|")
    if len(f) < 2:
        return None
    gid = f[-1]
    if gid.startswith(("RS_", "GB_")):
        gid = gid[3:]
    return f[0] + "__" + gid


want = set(cols)
keep_labels = []
for lf in t.leaf_nodes():
    lab = lf.taxon.label if lf.taxon else None
    if lab is not None and tip_key(lab) in want:
        keep_labels.append(lab)

print("tips matching the %d genomes: %d" % (len(cols), len(keep_labels)))
if len(keep_labels) != len(cols):
    found = {tip_key(l) for l in keep_labels}
    for g in sorted(want - found):
        print("  not in tree:", g)
    sys.exit("STOP: tip/genome mismatch.")

allb = [x.label for x in t.taxon_namespace if x.label is not None]
t.prune_taxa_with_labels([b for b in allb if b not in set(keep_labels)])
if hasattr(t, "purge_taxon_namespace"):
    t.purge_taxon_namespace()
t.suppress_unifurcations()
t.ladderize(ascending=False)
leaves = t.leaf_nodes()
print("subtree tips after prune: %d" % len(leaves))
if len(leaves) != len(cols):
    sys.exit("STOP: prune changed tip count.")

ycount = [0]
ymap = {}
for nd in t.postorder_node_iter():
    if nd.is_leaf():
        ymap[id(nd)] = float(ycount[0])
        ycount[0] += 1
    else:
        ch = [ymap[id(c)] for c in nd.child_nodes()]
        ymap[id(nd)] = (min(ch) + max(ch)) / 2.0
xmap = {}
for nd in t.preorder_node_iter():
    e = nd.edge.length
    e = 0.0 if e is None else float(e)
    xmap[id(nd)] = e if nd.parent_node is None else xmap[id(nd.parent_node)] + e

ordered = sorted(leaves, key=lambda z: ymap[id(z)])
row_order = [tip_key(lf.taxon.label) for lf in ordered]
ridx = [cols.index(g) for g in row_order]

grp_seq = [meta[g][1] for g in row_order]
runs = 1
for a, b in zip(grp_seq, grp_seq[1:]):
    if a != b:
        runs += 1
print("group runs down the tree: %d" % runs)


def order_columns(M):
    X = M[:, ridx]
    if X.shape[0] < 3:
        return np.arange(X.shape[0])
    D = pdist(X, metric="jaccard")
    return leaves_list(linkage(D, method="average"))


H1 = M1[order_columns(M1)][:, ridx].T
H2 = M2[order_columns(M2)][:, ridx].T
cmap = ListedColormap([ABSENT, PRESENT])

n = len(row_order)
total_units = len(units_og) + len(units_mm)
w1 = max(3.0, 9.0 * len(units_og) / max(total_units, 1))
w2 = max(3.0, 9.0 * len(units_mm) / max(total_units, 1))
fig = plt.figure(figsize=(4.0 + w1 + w2, max(8, 0.16 * n)))
gs = gridspec.GridSpec(
    1, 5, width_ratios=[2.2, 0.14, 0.14, w1, w2], wspace=0.03)

ax_t = fig.add_subplot(gs[0])
ax_c = fig.add_subplot(gs[1])
ax_g = fig.add_subplot(gs[2])
ax_1 = fig.add_subplot(gs[3])
ax_2 = fig.add_subplot(gs[4])

for nd in t.preorder_node_iter():
    y = ymap[id(nd)]
    x = xmap[id(nd)]
    if nd.parent_node is not None:
        ax_t.plot([xmap[id(nd.parent_node)], x], [y, y], color="k", lw=0.7)
    if not nd.is_leaf():
        ch = [ymap[id(c)] for c in nd.child_nodes()]
        ax_t.plot([x, x], [min(ch), max(ch)], color="k", lw=0.7)
ax_t.set_ylim(-0.5, n - 0.5)
ax_t.invert_yaxis()
ax_t.set_yticks([])
ax_t.set_xlabel("substitutions per site")
for s in ("top", "right", "left"):
    ax_t.spines[s].set_visible(False)

cat_colors = np.array([[CATC[meta[g][0]]] for g in row_order])
grp_colors = np.array([[GRPC[meta[g][1]]] for g in row_order])
ax_c.imshow(cat_colors, aspect="auto", interpolation="nearest")
ax_g.imshow(grp_colors, aspect="auto", interpolation="nearest")
for ax, lab in ((ax_c, "catalog"), (ax_g, "group")):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(lab, rotation=90, fontsize=8)

ax_1.imshow(H1, aspect="auto", interpolation="nearest", cmap=cmap,
            vmin=0, vmax=1)
ax_1.set_xlabel("%d eggNOG Bacteria-level OGs" % len(units_og))
ax_2.imshow(H2, aspect="auto", interpolation="nearest", cmap=cmap,
            vmin=0, vmax=1)
ax_2.set_xlabel("%d MMseqs protein clusters" % len(units_mm))
for ax in (ax_1, ax_2):
    ax.set_xticks([])
    ax.set_yticks([])

handles = [Patch(color=GRPC[g], label=GRPL[g])
           for g in ("amphibian", "non_amphibian", "context")]
handles += [Patch(color=CATC[c], label=CATL[c])
            for c in ("herptile", "ehi_amphibian", "youngblut", "gtdb_ref")]
handles += [Patch(color=PRESENT, label="present"),
            Patch(color=ABSENT, label="absent")]
fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
           fontsize=8)
fig.suptitle(
    "Angelakisella neighborhood, %d genomes: gene content in two "
    "currencies, rows ordered by the bac120 phylogeny" % n, fontsize=12)
fig.savefig(STEM + ".pdf", bbox_inches="tight")
fig.savefig(STEM + ".png", dpi=200, bbox_inches="tight")
print("wrote:", STEM + ".pdf")
print("wrote:", STEM + ".png")

print("")
print("CAPTION REQUIREMENTS this figure implies:")
print("- 61 genomes: 21 wild herptile, 10 EHI newt, 25 GTDB reference and")
print("  1 Youngblut Angelakisella, 4 Heteroruminococcus context")
print("- rows: bac120 subtree (MRCA of Angelakisella tips + 1 level)")
print("- left heatmap: eggNOG Bacteria-level OG presence")
print("- right heatmap: raw MMseqs cluster presence, same row order")
print("- columns Jaccard average-linkage within each panel")
print("- units in >= 3 of the 61 genomes")
print("- panel contrast = divergence splitting; cite the DIAMOND result")
print("  (66/67 amphibian-specific clusters had reference homologs)")

# FIG_ANGELAKISELLA_NEIGHBORHOOD_HEATMAP_V1
