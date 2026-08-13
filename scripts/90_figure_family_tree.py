# The family tree is drawn with tracks for source set, genus, relative
# evolutionary divergence, genus assignment and the assembly standard.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/fig_family_tree_red.py
# Output: results/figures/Figure_family_tree_red.pdf and .png
import os, sys, csv
from collections import Counter, OrderedDict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
import matplotlib.cm as cm
import dendropy

BASE = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
TREE = os.path.join(BASE, "work/rep_tree/figure_tree.nwk")
META = os.path.join(BASE, "work/rep_tree/figure_tree_metadata_genus.tsv")
WILD = os.path.join(BASE, "results/gtdbtk_wild_sgb_r220/gtdbtk.bac120.summary.tsv")
AMPH = os.path.join(BASE, "results/gtdbtk_ehi_amphibian_r220/gtdbtk.bac120.summary.tsv")
COH  = os.path.join(BASE, "results/unassigned_clade_coherence.tsv")
AUD  = os.path.join(BASE, "results/seqcode_representative_audit_v2.tsv")
FIGD = os.path.join(BASE, "results/figures")
TSV  = os.path.join(BASE, "results/family_tree_red_annotation.tsv")

def die(m):
    print("")
    print("!" * 72)
    print("FAILED: " + m)
    print("!" * 72)
    sys.exit(1)

def fnum(x):
    try:
        return float(str(x).strip())
    except Exception:
        return None

def norm(s):
    s = str(s).strip()
    if "|" in s:
        s = s.split("|")[-1]
    return s

print("=" * 72)
print("INPUTS")
print("=" * 72)
for p in (TREE, META, WILD):
    if not os.path.exists(p):
        die("missing " + p)
    print("  ok  " + p)

mrows = list(csv.DictReader(open(META), delimiter="\t"))
meta = {}
for r in mrows:
    meta[r["tip"].strip()] = {"arm": r["arm"].strip(),
                              "genus": (r["genus"] or "").strip(),
                              "genome": norm(r["genome"])}
print("  metadata rows: %d" % len(mrows))

red, gtdbgenus = {}, {}
for path, lab in ((WILD, "wild"), (AMPH, "newt")):
    if not os.path.exists(path):
        print("  ABSENT: %s" % path)
        continue
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    if "red_value" not in rows[0]:
        die("no red_value in " + path)
    n = 0
    for r in rows:
        g = (r["user_genome"] or "").strip()
        v = fnum(r["red_value"])
        if v is not None:
            red[g] = v
            n += 1
        gg = ""
        for f in (r["classification"] or "").split(";"):
            if f.strip().startswith("g__"):
                gg = f.strip()[3:].strip()
        gtdbgenus[g] = gg
    print("  %s: %d rows, %d with RED" % (lab, len(rows), n))

pass4 = set()
if os.path.exists(AUD):
    for r in csv.DictReader(open(AUD), delimiter="\t"):
        if (r.get("pass_four") or "0").strip() == "1":
            pass4.add((r["representative"] or "").strip())
print("  pass_four representatives: %d" % len(pass4))

deep = {}
if os.path.exists(COH):
    for r in csv.DictReader(open(COH), delimiter="\t"):
        c = (r["clade"] or "").strip()
        if c == "singleton":
            continue
        for g in (r["genomes"] or "").split(";"):
            if g.strip():
                deep[g.strip()] = c
print("  genomes in a multi-tip unassigned clade: %d" % len(deep))

tree = dendropy.Tree.get(path=TREE, schema="newick", preserve_underscores=True)
tree.ladderize(ascending=True)
tree.calc_node_root_distances(return_leaf_distances_only=False)
leaves = list(tree.leaf_node_iter())
labels = [nd.taxon.label if nd.taxon is not None else "" for nd in leaves]
n = len(leaves)
print("  tips: %d" % n)
if sum(1 for L in labels if L in meta) != n:
    die("metadata does not cover every tip")

QARM = ("herptile", "ehi_amphibian")
qtips = [L for L in labels if meta[L]["arm"] in QARM]
withred = [L for L in qtips if meta[L]["genome"] in red]
print("  query tips: %d   with a RED value: %d (%.1f%%)"
      % (len(qtips), len(withred), 100.0 * len(withred) / max(len(qtips), 1)))
if len(withred) < 0.8 * len(qtips):
    ex = [meta[L]["genome"] for L in qtips if meta[L]["genome"] not in red][:5]
    print("  unmatched examples: %s" % ex)
    print("  summary keys      : %s" % list(red)[:5])
    die("RED join covers under 80 percent of query tips")

vals = [red[meta[L]["genome"]] for L in withred]
vmin, vmax = min(vals), max(vals)
print("  RED range on tree: %.4f to %.4f" % (vmin, vmax))
assigned = [red[meta[L]["genome"]] for L in withred
            if gtdbgenus.get(meta[L]["genome"], "")]
unass = [red[meta[L]["genome"]] for L in withred
         if not gtdbgenus.get(meta[L]["genome"], "")]
print("  query tips WITH a GTDB genus   : %d" % len(assigned))
print("  query tips WITHOUT a GTDB genus: %d" % len(unass))
if assigned and unass:
    print("  median RED assigned %.4f   unassigned %.4f"
          % (float(np.median(assigned)), float(np.median(unass))))

ypos = {}
for i, nd in enumerate(leaves):
    ypos[nd] = float(i)
for nd in tree.postorder_node_iter():
    if not nd.is_leaf():
        ch = [ypos[c] for c in nd.child_node_iter()]
        ypos[nd] = sum(ch) / float(len(ch))
segs = []
for nd in tree.preorder_node_iter():
    p = nd.parent_node
    if p is not None:
        x0 = p.root_distance or 0.0
        segs.append([(x0, ypos[nd]), (nd.root_distance or x0, ypos[nd])])
for nd in tree.preorder_node_iter():
    kids = list(nd.child_node_iter())
    if kids:
        x = nd.root_distance or 0.0
        ys = [ypos[c] for c in kids]
        segs.append([(x, min(ys)), (x, max(ys))])
xmax = max((nd.root_distance or 0.0) for nd in tree.preorder_node_iter())

ARMCOL = OrderedDict([("herptile", "#b2182b"), ("ehi_amphibian", "#e08214"),
                      ("ehi", "#2166ac"), ("youngblut", "#762a83"),
                      ("gtdb_ref", "#c8d3de"), ("outgroup", "#000000")])
ARMLAB = {"herptile": "UHM herptile", "ehi_amphibian": "EHI newt", "ehi": "EHI mammal",
          "youngblut": "Youngblut", "gtdb_ref": "GTDB reference", "outgroup": "outgroup"}
seen = Counter(meta[L]["arm"] for L in labels)
gcount = Counter(meta[L]["genus"] for L in labels if meta[L]["genus"])
topg = [g for g, _ in gcount.most_common(14)]
PAL = ["#1f78b4", "#a6cee3", "#ff7f00", "#fdbf6f", "#33a02c", "#b2df8a", "#6a3d9a",
       "#cab2d6", "#8c510a", "#d8b365", "#01665e", "#80cdc1", "#666666", "#bf812d"]
GCOL = {g: PAL[i % len(PAL)] for i, g in enumerate(topg)}
NA = "#f2f2f2"
cmap = cm.get_cmap("viridis")

def c_arm(L):
    return ARMCOL.get(meta[L]["arm"], "#999999")

def c_genus(L):
    if meta[L]["arm"] == "outgroup":
        return "#ffffff"
    g = meta[L]["genus"]
    return "#111111" if not g else GCOL.get(g, "#d9d9d9")

def c_red(L):
    v = red.get(meta[L]["genome"])
    if v is None:
        return NA
    return cmap(1.0 - (v - vmin) / (vmax - vmin))

def c_unassigned(L):
    g = meta[L]["genome"]
    if meta[L]["arm"] not in QARM:
        return NA
    if g in deep:
        return "#000000"
    return "#bdbdbd" if not gtdbgenus.get(g, "") else "#ffffff"

def c_pass4(L):
    if meta[L]["arm"] not in QARM:
        return NA
    return "#1b7837" if meta[L]["genome"] in pass4 else "#ffffff"

fig = plt.figure(figsize=(16.0, max(14.0, n / 78.0)))
TL, TW, CW, GP = 0.050, 0.41, 0.016, 0.004
ax = fig.add_axes([TL, 0.045, TW, 0.925])
ax.add_collection(LineCollection(segs, colors="#3a3a3a", linewidths=0.30))
ax.set_xlim(-0.01 * xmax, xmax * 1.02)
ax.set_ylim(-1, n)
ax.invert_yaxis()
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.set_yticks([])
ax.tick_params(axis="x", labelsize=9)
ax.set_xlabel("substitutions per site\n(bac120, FastTree LG, no gamma)", fontsize=10)

cols = [("arm", c_arm), ("genus", c_genus),
        ("RED (dark = deeper)", c_red),
        ("no GTDB genus", c_unassigned),
        ("meets assembly standard", c_pass4)]
x = TL + TW + 0.012
for name, fn in cols:
    a2 = fig.add_axes([x, 0.045, CW, 0.925])
    for i in range(n):
        a2.add_patch(Rectangle((0, i), 1, 1, facecolor=fn(labels[i]), edgecolor="none"))
    a2.set_xlim(0, 1)
    a2.set_ylim(-1, n)
    a2.invert_yaxis()
    a2.set_xticks([])
    a2.set_yticks([])
    for s in ("top", "right", "left", "bottom"):
        a2.spines[s].set_visible(False)
    a2.set_xlabel(name, fontsize=7.5, rotation=90, ha="center", va="top", labelpad=6)
    x += CW + GP

runs, cur, st = [], None, 0
for i, L in enumerate(labels):
    g = meta[L]["genus"]
    if g != cur:
        if cur:
            runs.append((cur, st, i - 1))
        cur, st = g, i
if cur:
    runs.append((cur, st, n - 1))
axl = fig.add_axes([x, 0.045, 0.001, 0.925])
axl.set_ylim(-1, n)
axl.invert_yaxis()
axl.axis("off")
for g, a, b in runs:
    if (b - a + 1) >= 9:
        nh = sum(1 for i in range(a, b + 1) if meta[labels[i]]["arm"] == "herptile")
        lab = "%s (%d)" % (g, b - a + 1) + (", %d herptile" % nh if nh else "")
        axl.text(1.0, (a + b) / 2.0, lab, fontsize=7.5, va="center", ha="left",
                 color=GCOL.get(g, "#555555"))

cax = fig.add_axes([0.795, 0.055, 0.012, 0.16])
grad = np.linspace(1, 0, 256).reshape(-1, 1)
cax.imshow(grad, aspect="auto", cmap=cmap, extent=[0, 1, vmin, vmax])
cax.set_xticks([])
cax.yaxis.tick_right()
cax.tick_params(labelsize=7.5)
cax.set_title("RED", fontsize=8.5, pad=5)

h = [Patch(facecolor=ARMCOL[k], label="%s (%d)" % (ARMLAB[k], seen[k]))
     for k in ARMCOL if seen.get(k, 0) > 0]
h += [Patch(facecolor="#000000", label="in a genus-unassigned clade (%d)" % len(deep)),
      Patch(facecolor="#bdbdbd", label="no GTDB genus, ungrouped"),
      Patch(facecolor="#1b7837", label="meets 4 assembly recommendations (%d)" % len(pass4)),
      Patch(facecolor=NA, label="not applicable (reference / outgroup)")]
fig.legend(handles=h, loc="upper right", bbox_to_anchor=(0.995, 0.96),
           frameon=False, fontsize=9)
fig.suptitle("Ruminococcaceae, %d genomes across five catalogs, with relative "
             "evolutionary divergence of query genomes" % n, fontsize=13.5, y=0.985)

for e in ("pdf", "png"):
    o = os.path.join(FIGD, "Figure_family_tree_red." + e)
    fig.savefig(o, dpi=300, bbox_inches="tight")
    print("  wrote " + o)
plt.close(fig)

if os.path.exists(TSV):
    print("  NOT overwriting existing " + TSV)
else:
    with open(TSV, "w") as fh:
        fh.write("tip\tarm\tgenus_on_tree\tgenome\tgtdb_genus\tred\tunassigned_clade\tpass_four\ty_order\n")
        for i, L in enumerate(labels):
            g = meta[L]["genome"]
            fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\n"
                     % (L, meta[L]["arm"], meta[L]["genus"], g,
                        gtdbgenus.get(g, ""), red.get(g, ""),
                        deep.get(g, ""), int(g in pass4), i))
    print("  wrote " + TSV)

print("")
print("=" * 72)
print("CAPTION REQUIREMENTS")
print("=" * 72)
print("  1. RED IS ONLY DEFINED FOR QUERY GENOMES. It is computed during")
print("     GTDB-Tk placement on the GTDB reference tree, NOT on this tree.")
print("     Reference and outgroup tips are 'not applicable', never low.")
print("  2. RED IS A CONTINUOUS RAMP HERE ON PURPOSE. Do not describe any tip")
print("     as passing or failing a RED threshold. A 5th-percentile cut split")
print("     an 11-tip clade whose internal RED range was 0.007.")
print("  3. RED and the branch lengths are INDEPENDENT. The ramp is not a")
print("     recolouring of the x axis.")
print("  4. The unassigned-clade column is GTDB-Tk assignment plus monophyly.")
print("     NO patristic cutoff was used anywhere in this figure.")
print("  5. Reference arm is the UNFILTERED 1,247. Tips are GENOMES, not the")
print("     dRep clusters used as units in the composition results.")
print("  6. EHI mammal is Rodentia/Carnivora/Lagomorpha plus 2 marsupials and")
print("     1 parrot. Never call it endotherm gut.")
print("  7. Absence of herptile tips from a clade is NON-RECOVERY, not absence")
print("     from herptile guts.")
print("")
print("FIG_FAMILY_TREE_RED_V1_20260806 COMPLETE")
# FIG_FAMILY_TREE_RED_V1_20260806
