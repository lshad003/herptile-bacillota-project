# The CAZy family heatmap is drawn against a pruned subtree of the focal genomes.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/fig_cazy_tree_heatmap_v2.py
# Output: results/figures/Figure_cazy_tree_heatmap_v2.pdf and .png,
#         results/cazy_tree_heatmap_v2_order.tsv
import os, sys, csv
from collections import Counter, OrderedDict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle, Patch
import dendropy

BASE = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
TREE = os.path.join(BASE, "work/rep_tree/figure_tree.nwk")
META = os.path.join(BASE, "work/rep_tree/figure_tree_metadata_genus.tsv")
MTX  = os.path.join(BASE, "results/cazy_focal_family_matrix.tsv")
HAP  = os.path.join(BASE, "results/happi_cazy_focal.tsv")
FIGD = os.path.join(BASE, "results/figures")
OUTT = os.path.join(BASE, "results/cazy_tree_heatmap_v2_order.tsv")
Q = 0.05

def die(m):
    print("")
    print("!" * 72)
    print("FAILED: " + m)
    print("!" * 72)
    sys.exit(1)

def norm(s):
    s = str(s).strip()
    return s.split("|")[-1] if "|" in s else s

def cls(f):
    for p in ("GH", "GT", "PL", "CE", "AA", "CBM"):
        if f.startswith(p):
            return p
    return "other"

print("=" * 72)
print("STEP 1  INPUTS")
print("=" * 72)
for p in (TREE, META, MTX, HAP):
    if not os.path.exists(p):
        die("missing " + p)
rows = list(csv.DictReader(open(MTX), delimiter="\t"))
cz = {r["genome_id"]: r for r in rows}
meta_cols = {"genome_id", "set", "genus", "arm", "clade", "total"}
allfams = [c for c in rows[0].keys() if c not in meta_cols]
print("  genomes: %d   families detected: %d" % (len(cz), len(allfams)))
print("  by set: %s" % dict(Counter(r["set"] for r in rows)))

hp = list(csv.DictReader(open(HAP), delimiter="\t"))
tested = {r["family"] for r in hp}
sig, direction = {}, {}
for r in hp:
    if r["q"] not in ("", "NA") and float(r["q"]) < Q:
        sig[r["family"]] = float(r["q"])
        direction[r["family"]] = 1 if float(r["diff"]) > 0 else -1
print("  tested: %d   significant at q<%.2f: %d" % (len(tested), Q, len(sig)))
print("    amphibian-higher %d   reference-higher %d"
      % (sum(1 for v in direction.values() if v > 0),
         sum(1 for v in direction.values() if v < 0)))

focal = [g for g in cz if cz[g]["set"] == "focal125"]
unass = [g for g in cz if cz[g]["set"] == "unassigned_clade"]
print("  focal125: %d   unassigned_clade: %d" % (len(focal), len(unass)))

print("")
print("=" * 72)
print("STEP 2  COLUMN ORDER")
print("=" * 72)
prev_all = {f: sum(1 for g in cz if int(cz[g][f]) > 0) for f in allfams}
def rank(f):
    s = 0 if f in sig and direction[f] > 0 else (2 if f in sig else 1)
    return (s, cls(f), -prev_all[f], f)
fams = sorted(allfams, key=rank)
print("  columns ordered: amphibian-higher, then untested/ns, then reference-higher")
print("  by class: %s" % dict(Counter(cls(f) for f in fams).most_common()))
print("  NOT tested (too rare or near-universal): %d" % len(set(allfams) - tested))

print("")
print("=" * 72)
print("STEP 3  PRUNE TO THE FOCAL SET ONLY")
print("=" * 72)
mrows = list(csv.DictReader(open(META), delimiter="\t"))
tipof = {norm(r["genome"]): r["tip"].strip() for r in mrows}
want = {g: tipof[g] for g in focal if g in tipof}
if len(want) != len(focal):
    die("not every focal genome has a tip")
tree = dendropy.Tree.get(path=TREE, schema="newick", preserve_underscores=True)
tree.retain_taxa_with_labels(list(want.values()))
tree.ladderize(ascending=True)
tree.calc_node_root_distances(return_leaf_distances_only=False)
leaves = list(tree.leaf_node_iter())
labels = [nd.taxon.label for nd in leaves if nd.taxon is not None]
n = len(labels)
print("  pruned tips: %d" % n)
print("  THE UNASSIGNED CLADES ARE NOT ON THIS TREE. In v1 they were, and")
print("  pruning collapsed 8 clades with different sisters into one apparent")
print("  block. They are drawn as a separate panel below instead.")
tip2g = {v: k for k, v in want.items()}

M = np.zeros((n, len(fams)), dtype=np.int8)
for i, L in enumerate(labels):
    r = cz[tip2g[L]]
    for j, f in enumerate(fams):
        M[i, j] = 1 if int(r[f]) > 0 else 0
print("  focal matrix: %d x %d   density %.3f" % (M.shape[0], M.shape[1], M.mean()))

unass.sort(key=lambda g: (cz[g]["clade"] == "singleton", cz[g]["clade"], g))
U = np.zeros((len(unass), len(fams)), dtype=np.int8)
for i, g in enumerate(unass):
    for j, f in enumerate(fams):
        U[i, j] = 1 if int(cz[g][f]) > 0 else 0
print("  unassigned matrix: %d x %d   density %.3f" % (U.shape[0], U.shape[1], U.mean()))

ypos = {}
for i, nd in enumerate(leaves):
    ypos[nd] = float(i)
for nd in tree.postorder_node_iter():
    if not nd.is_leaf():
        k = list(nd.child_node_iter())
        ypos[nd] = sum(ypos[c] for c in k) / float(len(k))
segs = []
for nd in tree.preorder_node_iter():
    p = nd.parent_node
    if p is not None:
        x0 = p.root_distance or 0.0
        segs.append([(x0, ypos[nd]), (nd.root_distance or x0, ypos[nd])])
for nd in tree.preorder_node_iter():
    k = list(nd.child_node_iter())
    if k:
        x = nd.root_distance or 0.0
        ys = [ypos[c] for c in k]
        segs.append([(x, min(ys)), (x, max(ys))])
xmax = max((nd.root_distance or 0.0) for nd in tree.preorder_node_iter())

ARM = {"amphibian": "#b2182b", "reference": "#2166ac", "endotherm": "#000000",
       "herptile": "#b2182b", "ehi_amphibian": "#e08214"}
GEN = {"Anaerotruncus": "#1f78b4", "UBA866": "#a6cee3"}
CLADE = ["#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
         "#b3de69", "#fccde5"]

# Panel heights are proportional to row count so that a row in the focal
# block and a row in the unassigned block are drawn at the same height.
N_TOP, N_BOT = M.shape[0], U.shape[0]
ROW_IN = 0.105
GAP_IN = 0.85
MARGIN_TOP_IN, MARGIN_BOT_IN = 1.1, 1.6
FIG_W = 19.0
FIG_H = (MARGIN_BOT_IN + N_BOT * ROW_IN + GAP_IN
         + N_TOP * ROW_IN + MARGIN_TOP_IN)

fig = plt.figure(figsize=(FIG_W, FIG_H))
H_TOP_FRAC = (N_TOP * ROW_IN) / FIG_H
H_BOT_FRAC = (N_BOT * ROW_IN) / FIG_H
HBOT = MARGIN_BOT_IN / FIG_H
HTOP = HBOT + H_BOT_FRAC + (GAP_IN / FIG_H)
TL, TW, CW, GP = 0.040, 0.155, 0.011, 0.003
HX, HW = TL + TW + 0.010 + 2 * (CW + GP), 0.60
print("layout: %d top rows, %d bottom rows, figure %.1f x %.1f inches"
      % (N_TOP, N_BOT, FIG_W, FIG_H))

ax = fig.add_axes([TL, HTOP, TW, H_TOP_FRAC])
ax.add_collection(LineCollection(segs, colors="#3a3a3a", linewidths=0.45))
ax.set_xlim(-0.01 * xmax, xmax * 1.02)
ax.set_ylim(-1, n)
ax.invert_yaxis()
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.set_yticks([])
ax.tick_params(axis="x", labelsize=8)
ax.set_xlabel("substitutions per site (bac120 subtree)", fontsize=8)

x = TL + TW + 0.010
for name, fn in [("arm", lambda L: ARM.get(cz[tip2g[L]]["arm"], "#999")),
                 ("genus", lambda L: GEN.get(cz[tip2g[L]]["genus"], "#e0e0e0"))]:
    a2 = fig.add_axes([x, HTOP, CW, H_TOP_FRAC])
    for i, L in enumerate(labels):
        a2.add_patch(Rectangle((0, i), 1, 1, facecolor=fn(L), edgecolor="none"))
    a2.set_xlim(0, 1)
    a2.set_ylim(-1, n)
    a2.invert_yaxis()
    a2.set_xticks([])
    a2.set_yticks([])
    for s in ("top", "right", "left", "bottom"):
        a2.spines[s].set_visible(False)
    a2.text(0.5, -0.004, name, transform=a2.transAxes, fontsize=7.5,
            rotation=90, ha="center", va="top")
    x += CW + GP

cmap = matplotlib.colors.ListedColormap(["#ececec", "#1f4e79"])
axh = fig.add_axes([HX, HTOP, HW, H_TOP_FRAC])
axh.imshow(M, aspect="auto", interpolation="nearest", vmin=0, vmax=1, cmap=cmap)
STEP_A = max(1, len(fams) // 90)
axh.set_xticks(range(0, len(fams), STEP_A))
axh.set_xticklabels([fams[k] for k in range(0, len(fams), STEP_A)],
                    rotation=90, fontsize=5.5)
axh.set_yticks([])
axh.set_title("a   %d focal genomes, Anaerotruncus and UBA866. The %d tested are the amphibian and GTDB reference genomes; the endotherm-labelled genome is drawn but excluded from the test." % (M.shape[0], M.shape[0] - 1),
              loc="left", fontsize=11, fontweight="bold")
for s in ("top", "right", "left", "bottom"):
    axh.spines[s].set_visible(False)

axm = fig.add_axes([HX, HTOP - (0.34 / FIG_H), HW, 0.24 / FIG_H])
for j, f in enumerate(fams):
    c = "#b2182b" if f in sig and direction[f] > 0 else \
        ("#2166ac" if f in sig else ("#bdbdbd" if f in tested else "#ffffff"))
    axm.add_patch(Rectangle((j, 0), 1, 1, facecolor=c, edgecolor="none"))
axm.set_xlim(0, len(fams))
axm.set_ylim(0, 1)
axm.set_xticks([])
axm.set_yticks([])
for s in ("top", "right", "left", "bottom"):
    axm.spines[s].set_visible(False)
axm.text(-0.004, 0.5, "test", transform=axm.transAxes, fontsize=7.5,
         rotation=0, ha="right", va="center")

axu = fig.add_axes([HX, HBOT, HW, H_BOT_FRAC])
axu.imshow(U, aspect="auto", interpolation="nearest", vmin=0, vmax=1, cmap=cmap)
axu.set_yticks([])
step = max(1, len(fams) // 90)
axu.set_xticks(range(0, len(fams), step))
axu.set_xticklabels([fams[k] for k in range(0, len(fams), step)],
                    rotation=90, fontsize=5.5)
for s in ("top", "right", "left", "bottom"):
    axu.spines[s].set_visible(False)
axu.set_title("b   29 genus-unassigned genomes, DESCRIPTIVE ONLY, not tested. "
              "Row order is by clade, NOT phylogeny: these 8 clades sit in "
              "different parts of the family tree with different sister genera.",
              loc="left", fontsize=10, fontweight="bold")
axu.set_xlabel("%d CAZy families detected across all 154 genomes, ordered: "
               "amphibian-higher, untested or not significant, reference-higher"
               % len(fams), fontsize=9)

axc = fig.add_axes([HX - (CW + GP), HBOT, CW, H_BOT_FRAC])
for i, g in enumerate(unass):
    c = cz[g]["clade"]
    if c == "singleton":
        col = "#000000"
    else:
        try:
            col = CLADE[(int(c) - 1) % len(CLADE)]
        except Exception:
            col = "#ffffff"
    axc.add_patch(Rectangle((0, i), 1, 1, facecolor=col, edgecolor="none"))
axc.set_xlim(0, 1)
axc.set_ylim(0, len(unass))
axc.invert_yaxis()
axc.set_xticks([])
axc.set_yticks([])
for s in ("top", "right", "left", "bottom"):
    axc.spines[s].set_visible(False)
axc.text(0.5, -0.004, "clade", transform=axc.transAxes, fontsize=7.5,
         rotation=90, ha="center", va="top")

h = [Patch(facecolor=ARM["amphibian"], label="amphibian / herptile"),
     Patch(facecolor=ARM["reference"], label="GTDB reference"),
     Patch(facecolor=ARM["ehi_amphibian"], label="EHI newt"),
     Patch(facecolor=ARM["endotherm"], label="endotherm-labelled, excluded from test"),
     Patch(facecolor=GEN["Anaerotruncus"], label="Anaerotruncus"),
     Patch(facecolor=GEN["UBA866"], label="UBA866 / Paludihabitans"),
     Patch(facecolor="#1f4e79", label="family present"),
     Patch(facecolor="#ececec", label="family absent"),
     Patch(facecolor="#b2182b", label="q<0.05, higher in amphibian (%d)"
           % sum(1 for v in direction.values() if v > 0)),
     Patch(facecolor="#2166ac", label="q<0.05, higher in reference (%d)"
           % sum(1 for v in direction.values() if v < 0)),
     Patch(facecolor="#bdbdbd", label="tested, not significant (%d)"
           % (len(tested) - len(sig))),
     Patch(facecolor="#ffffff", edgecolor="#999999",
           label="not tested, too rare or near-universal (%d)"
           % (len(fams) - len(tested)))]
fig.legend(handles=h, loc="upper left",
           bbox_to_anchor=(0.035, min(1.0, HTOP + H_TOP_FRAC + (1.15 / FIG_H))),
           frameon=False, fontsize=8.5, ncol=4)
fig.suptitle("CAZy family presence: within-genus amphibian versus reference contrast, "
             "all %d detected families" % len(fams), fontsize=13,
             y=min(0.995, HTOP + H_TOP_FRAC + (0.62 / FIG_H)))

for e in ("pdf", "png"):
    o = os.path.join(FIGD, "Figure_cazy_tree_heatmap_v2." + e)
    fig.savefig(o, dpi=250, bbox_inches="tight")
    print("  wrote " + o)
plt.close(fig)

if os.path.exists(OUTT):
    print("  NOT overwriting existing " + OUTT)
else:
    with open(OUTT, "w") as fh:
        fh.write("panel\torder\tgenome_id\tset\tarm\tgenus\tclade\tn_families\n")
        for i, L in enumerate(labels):
            r = cz[tip2g[L]]
            fh.write("a\t%d\t%s\t%s\t%s\t%s\t%s\t%d\n"
                     % (i, tip2g[L], r["set"], r["arm"], r["genus"], r["clade"],
                        int(M[i].sum())))
        for i, g in enumerate(unass):
            r = cz[g]
            fh.write("b\t%d\t%s\t%s\t%s\t%s\t%s\t%d\n"
                     % (i, g, r["set"], r["arm"], r["genus"], r["clade"],
                        int(U[i].sum())))
    print("  wrote " + OUTT)

print("")
print("=" * 72)
print("CAPTION REQUIREMENTS")
print("=" * 72)
print("  1. THE TEST IS WITHIN GENUS. Model is presence ~ amphibian + genus,")
print("     so genus is adjusted for and the amphibian term is the contrast:")
print("     46 vs 16 in Anaerotruncus, 52 vs 10 in UBA866.")
print("  2. Only these two genera have amphibian and reference genomes")
print("     phylogenetically interleaved (Slatkin-Maddison, >=5 transitions).")
print("     Elsewhere host and lineage are the same variable.")
print("  3. Panel a is a PRUNED SUBTREE, not a new inference.")
print("  4. Panel b rows are NOT phylogenetically ordered and the clades are")
print("     NOT sister to the focal genera. Sisters are Hydrogenoanaerobacterium,")
print("     Faecivivens, Harryflintia, Avimicrobium and UBA1405.")
print("  5. All %d detected families are shown. %d were tested; the rest were" % (len(fams), len(tested)))
print("     too rare or near-universal to fit.")
print("  6. dbCAN v13.0, i-Evalue < 1e-15, coverage > 0.35. Nulls 101 and 202")
print("     returned 0 of 196. 8 of 196 raised a spline boundary warning.")
print("  7. Presence/absence, not copy number.")
print("")
print("FIG_CAZY_TREE_HEATMAP_V2_20260806 COMPLETE")
# FIG_CAZY_TREE_HEATMAP_V2_20260806
