#!/usr/bin/env python3
# Figure 2 is drawn and the genus-expansion permutation is run at 9,999 iterations, seed 20260803.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/fig2_novelty_v6.py
# Output: results/figures/Figure2_novelty_expansion.pdf/.png, results/genus_expansion_permutation.tsv
import os, sys, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
NOV = B + "/results/novelty_proportions.tsv"
EXP = B + "/results/genus_expansion.tsv"
SGB = B + "/data/sgb_manifest.tsv"
TAX = "/srv/projects/db/gtdbtk/220/taxonomy/bac120_taxonomy_r220_reps.tsv"
OUTD = B + "/results/figures"
STEM = OUTD + "/Figure2_novelty_expansion"
PERM_OUT = B + "/results/genus_expansion_permutation.tsv"

FAM = "Ruminococcaceae"
NPERM = 9999
SEED = 20260803
UNK = ("UNASSIGNED", "(unassigned)", "")

for p in (NOV, EXP, SGB, TAX):
    if not os.path.exists(p):
        print("MISSING:", p); sys.exit(1)
if not os.path.isdir(OUTD):
    os.makedirs(OUTD)

nov = list(csv.DictReader(open(NOV), delimiter="\t"))
exp = list(csv.DictReader(open(EXP), delimiter="\t"))

hosts = {}
with open(SGB) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    I = {k: i for i, k in enumerate(h)}
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(I.values()):
            continue
        if p[I["family"]] != FAM or p[I["has_wild"]] != "yes":
            continue
        gname = p[I["genus"]]
        if gname in UNK:
            continue
        hosts.setdefault(gname, set()).update(
            x for x in p[I["host_species"]].split(";") if x)
print("genera with host data: %d" % len(hosts))

fig = plt.figure(figsize=(18.5, 6.8))
gs = fig.add_gridspec(1, 3, width_ratios=[0.95, 1.55, 0.90], wspace=0.62)

a = fig.add_subplot(gs[0])
short = ["all Bacillota_A\nall animals", "all Bacillota_A\nwild only",
         "Ruminococcaceae\nall animals", "Ruminococcaceae\nwild only"]
n = np.array([int(r["n_sgbs"]) for r in nov], dtype=float)
ns = np.array([int(r["novel_species"]) for r in nov], dtype=float)
ng = np.array([int(r["novel_genus"]) for r in nov], dtype=float)
nf = np.array([int(r["novel_family"]) for r in nov], dtype=float)
matched = n - ns
x = np.arange(len(nov))
bot = np.zeros(len(nov))
for v, colr, lab in ((matched, "#9E9E9E", "matched a GTDB species cluster"),
                     (ns - ng, "#7BA7D7", "novel species, known genus"),
                     (ng - nf, "#3E6DA8", "novel genus, known family"),
                     (nf, "#16294A", "novel family")):
    pc = 100.0 * v / n
    a.bar(x, pc, bottom=bot, width=0.62, color=colr, edgecolor="black",
          lw=0.5, label=lab)
    bot += pc
# v6: the "N matched" callouts were dropped BELOW the axis at -26 points,
# where they ran straight through the rotated tick labels. Put the count
# INSIDE the grey band instead, which is where the quantity actually lives.
for i in range(len(nov)):
    if matched[i] > 0:
        pcm = 100.0 * matched[i] / n[i]
        a.text(i, pcm / 2.0, "%d" % int(matched[i]), ha="center", va="center",
               fontsize=7.5, color="black")
    a.text(i, 101.5, "n=%d" % int(n[i]), ha="center", fontsize=8.5)
a.set_xticks(x)
a.set_xticklabels(short, fontsize=8.0, rotation=22, ha="right",
                  rotation_mode="anchor")
a.set_ylabel("percent of SGBs")
a.set_ylim(0, 108)
a.set_xlim(-0.6, len(nov) - 0.4)
a.set_title("a   Novelty against GTDB r220", loc="left", fontweight="bold",
            fontsize=12.5)
a.legend(fontsize=7.5, frameon=False, loc="lower left",
         bbox_to_anchor=(-0.02, 1.10), ncol=2)
a.grid(axis="y", alpha=0.2, lw=0.5)

b = fig.add_subplot(gs[1])
g = [r["genus"] for r in exp]
sx = np.array([int(r["gtdb_species_clusters"]) for r in exp], dtype=float)
sy = np.array([int(r["wild_sgbs"]) for r in exp], dtype=float)
nh = np.array([len(hosts.get(k, ())) for k in g], dtype=float)
rat = sy / np.maximum(sx, 1)
col = ["#C44E52" if v >= 2 else ("#4C72B0" if v <= 0.5 else "#9E9E9E")
       for v in rat]
size = np.clip(nh * 26, 30, 330)

lo, hi = 0.62, max(sx.max(), sy.max()) * 2.2
b.set_xscale("log"); b.set_yscale("log")
b.set_xlim(lo, hi); b.set_ylim(lo, hi)
b.plot([lo, hi], [lo, hi], color="black", ls="--", lw=1.2, zorder=1)
b.plot([lo, hi / 2.0], [2 * lo, hi], color="#C44E52", ls=":", lw=1.1, zorder=1)
b.scatter(sx, sy, s=size, c=col, edgecolor="black", lw=0.6, alpha=0.9, zorder=3)

show = [i for i in range(len(g))
        if rat[i] >= 2 or rat[i] <= 0.35 or sy[i] >= 15]
above = sorted([i for i in show if rat[i] >= 1.0], key=lambda i: -np.log10(sy[i]))
below = sorted([i for i in show if rat[i] < 1.0], key=lambda i: -np.log10(sy[i]))

def park(idxs, xfrac, ha):
    if not idxs:
        return
    top, botf = 0.97, 0.05
    step = (top - botf) / max(len(idxs) - 1, 1)
    for k, i in enumerate(idxs):
        yf = top - k * step
        b.annotate(g[i], xy=(sx[i], sy[i]), xycoords="data",
                   xytext=(xfrac, yf), textcoords="axes fraction",
                   fontsize=7.8, ha=ha, va="center",
                   annotation_clip=False,
                   arrowprops=dict(arrowstyle="-", color="#AAAAAA", lw=0.6,
                                   shrinkA=2, shrinkB=4))

park(above, -0.045, "right")
park(below, 1.045, "left")

H = np.log10(hi)
b.text(10 ** (H - 0.055), 10 ** (H - 0.02), "parity", fontsize=8.5,
       color="black", ha="right", va="top",
       bbox=dict(fc="white", ec="none", alpha=0.9, pad=0.8))
b.text(10 ** (H - 0.355), 10 ** (H - 0.02), "2x", fontsize=8.5,
       color="#C44E52", ha="right", va="top",
       bbox=dict(fc="white", ec="none", alpha=0.9, pad=0.8))

b.set_xlabel("GTDB r220 species clusters in genus")
# v6: labelpad was 118, which pushed this label out of panel b and into
# panel a as free-floating rotated text. The label column is parked in
# axes-fraction coords and does not occupy the y-label slot.
b.set_ylabel("wild-host SGBs recovered", labelpad=6)
b.set_title("b   Expansion is genus-specific", loc="left", fontweight="bold",
            fontsize=12.5)
b.grid(alpha=0.18, lw=0.5, which="both")

# v6: panel b had no colour legend at all. The top legend is panel a's.
n_red = int((rat >= 2).sum())
n_blue = int((rat <= 0.5).sum())
n_grey = len(rat) - n_red - n_blue
handles = [
    Line2D([], [], marker="o", ls="", mfc="#C44E52", mec="black", mew=0.6,
           ms=7, label="$\\geq$2x more SGBs (n=%d)" % n_red),
    Line2D([], [], marker="o", ls="", mfc="#9E9E9E", mec="black", mew=0.6,
           ms=7, label="between 0.5x and 2x (n=%d)" % n_grey),
    Line2D([], [], marker="o", ls="", mfc="#4C72B0", mec="black", mew=0.6,
           ms=7, label="$\\leq$0.5x (n=%d)" % n_blue),
]
b.legend(handles=handles, fontsize=7.5, frameon=True, framealpha=0.9,
         loc="lower right", bbox_to_anchor=(0.985, 0.015), handletextpad=0.4,
         borderpad=0.5)
b.text(0.5, -0.155, "point area: host species per genus (1 to 10)",
       transform=b.transAxes, ha="center", va="top", fontsize=7.8,
       color="#444444")

c = fig.add_subplot(gs[2])
def fld(t, pre):
    for z in t.split(";"):
        z = z.strip()
        if z.startswith(pre):
            return z[len(pre):]
    return ""
ref = {}
with open(TAX) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2 and fld(p[1], "f__") == FAM:
            gg = fld(p[1], "g__")
            if gg:
                ref[gg] = ref.get(gg, 0) + 1
gl = sorted(ref)
w = np.array([ref[k] for k in gl], dtype=float)
w = w / w.sum()
tot = int(sy.sum())
obs = int((rat >= 2.0).sum())
rng = np.random.RandomState(SEED)
null = np.zeros(NPERM, dtype=int)
for i in range(NPERM):
    dr = rng.multinomial(tot, w)
    null[i] = int(sum(1 for j in range(len(gl)) if dr[j] / ref[gl[j]] >= 2.0))
pv = (1 + int((null >= obs).sum())) / float(NPERM + 1)

# v6: write the test result to disk so the manuscript number can be checked
# against a file instead of being regenerated at plot time.
if os.path.exists(PERM_OUT):
    print("NOTE: %s exists, not rewriting it" % PERM_OUT)
else:
    with open(PERM_OUT, "w") as fh:
        fh.write("quantity\tvalue\n")
        fh.write("observed_genera_at_2x\t%d\n" % obs)
        fh.write("n_permutations\t%d\n" % NPERM)
        fh.write("seed\t%d\n" % SEED)
        fh.write("p_value\t%.6f\n" % pv)
        fh.write("wild_sgbs_allocated\t%d\n" % tot)
        fh.write("gtdb_genera_in_null\t%d\n" % len(gl))
        fh.write("null_source\t%s\n" % os.path.basename(TAX))
        fh.write("null_unit\tone row per GTDB r220 species representative\n")
        fh.write("null_max\t%d\n" % int(null.max()))
        fh.write("null_mean\t%.4f\n" % float(null.mean()))
    print("wrote %s" % PERM_OUT)

mx = max(int(null.max()), obs) + 1
c.hist(null, bins=np.arange(-0.5, mx + 1.5, 1), color="#C8C8C8",
       edgecolor="black", lw=0.5)
c.axvline(obs, color="#C44E52", lw=2.4, zorder=5)
yl = c.get_ylim()[1]
c.text(obs - 0.3, yl * 0.92, "observed = %d\np = %.4f" % (obs, pv),
       color="#C44E52", fontsize=10, fontweight="bold", ha="right", va="top")
c.set_xlabel("genera with $\\geq$2x more SGBs\nthan GTDB species clusters")
c.set_ylabel("permutations")
c.set_title("c   Database-proportional null", loc="left", fontweight="bold",
            fontsize=12.5)
c.grid(axis="y", alpha=0.2, lw=0.5)
c.text(0.97, 0.60, "%d wild SGBs allocated across\n%d GTDB genera in proportion\n"
       "to their species-cluster counts\n%d permutations"
       % (tot, len(gl), NPERM), transform=c.transAxes, ha="right", va="top",
       fontsize=7.5, color="#444444")

fig.savefig(STEM + ".png", dpi=300, bbox_inches="tight")
fig.savefig(STEM + ".pdf", bbox_inches="tight")
print("wrote", STEM + ".png")
print("labelled %d genera: %d left column, %d right column"
      % (len(show), len(above), len(below)))
print("panel b colours: %d red, %d grey, %d blue" % (n_red, n_grey, n_blue))
print("observed %d genera at >=2x, permutation p = %.4f" % (obs, pv))
print()
print("CAPTION MUST STATE:")
print("  point area encodes host species per genus (1 to 10)")
print("  panel b axes are log10; a genus at exactly 2.0x sits ON the dotted")
print("    line, not above it. Two of the %d red points are at exactly 2.0." % n_red)
print("  panel c null draws from %d GTDB r220 species representatives" % len(gl))
print("  the 2x threshold was chosen after inspecting the data")
# FIG2_NOVELTY_V6_20260808
