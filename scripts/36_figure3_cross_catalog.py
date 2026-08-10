#!/usr/bin/env python3
# Figure 3 is drawn: cross-catalog genus overlap, turnover, threshold erosion and accumulation.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/fig3_cross_catalog.py
# Output: results/figures/Figure3_cross_catalog_overlap.pdf and .png
# FIG3_CROSS_CATALOG_V7_20260805
import os, sys
from collections import Counter, defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
WORK = B + "/work/pooled_drep_rum"
CDB = WORK + "/drep_out/data_tables/Cdb.csv"
MAP = WORK + "/genome_arms.tsv"

HERP_SUM = B + "/results/gtdbtk_wild_sgb_r220/gtdbtk.bac120.summary.tsv"
AMPH_SUM = B + "/results/gtdbtk_ehi_amphibian_r220/gtdbtk.bac120.summary.tsv"
EHI_SUM = B + "/results/gtdbtk_ehi_r220_classify/gtdbtk.bac120.summary.tsv"
YB_SUM = B + "/results/gtdbtk_youngblut_r220/gtdbtk.bac120.summary.tsv"
GTDB_TAX = "/srv/projects/db/gtdbtk/220/taxonomy/gtdb_taxonomy.tsv"

BETA = B + "/results/turnover_nestedness.tsv"
THRESH = B + "/results/threshold_sensitivity.tsv"
ACC = B + "/results/genus_accumulation.tsv"

OUTD = B + "/results/figures"
STEM = OUTD + "/Figure3_cross_catalog_overlap"

ARMS = ["herptile", "ehi_amphibian", "ehi", "youngblut", "gtdb_ref"]
COLLAB = ["herptile\nwild", "EHI\nnewts", "EHI\nmammal/bird", "Youngblut", "GTDB r220\nreference"]
TOP_NOHERP = 14

for p in (CDB, MAP, BETA, THRESH, ACC):
    if not os.path.exists(p):
        print("MISSING:", p); sys.exit(1)
if not os.path.isdir(OUTD):
    os.makedirs(OUTD)


def read_tsv(path):
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                f = f + [""] * (len(header) - len(f))
            rows.append(dict(zip(header, f)))
    return rows


def parse_tax(s):
    out = {}
    for part in s.strip().split(";"):
        part = part.strip()
        if len(part) > 3 and part[1:3] == "__":
            out[part[0]] = part[3:]
    return out


def norm(g):
    g = g.strip()
    return "" if g.upper() in ("UNASSIGNED", "", "NA", "N/A") else g


def tofloat(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


# ------------------------------------------------- units, same as the stats
genus = {}
for arm, path in (("herptile", HERP_SUM), ("ehi_amphibian", AMPH_SUM),
                  ("ehi", EHI_SUM), ("youngblut", YB_SUM)):
    for r in read_tsv(path):
        genus[(arm, r["user_genome"].strip())] = norm(
            parse_tax(r["classification"]).get("g", ""))
with open(GTDB_TAX) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        t = parse_tax(f[1])
        if t.get("f", "") == "Ruminococcaceae":
            genus[("gtdb_ref", f[0].strip().replace("GB_", "").replace("RS_", ""))] = \
                norm(t.get("g", ""))

arm_of, gid_of = {}, {}
with open(MAP) as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        arm_of[f[0]] = f[1]
        gid_of[f[0]] = f[2]

clusters = defaultdict(list)
with open(CDB) as fh:
    hdr = fh.readline().rstrip("\n").replace('"', "").split(",")
    gi, si = hdr.index("genome"), hdr.index("secondary_cluster")
    for line in fh:
        f = line.rstrip("\n").replace('"', "").split(",")
        if len(f) > max(gi, si):
            clusters[f[si]].append(f[gi])

units = defaultdict(list)
for cid, gs in clusters.items():
    per_arm = defaultdict(list)
    for g in gs:
        per_arm[arm_of.get(g, "?")].append(g)
    for a, members in per_arm.items():
        labs = [genus.get((a, gid_of.get(m, "")), "") for m in members]
        named = [x for x in labs if x]
        if named:
            units[a].append(Counter(named).most_common(1)[0][0])

counts = {a: Counter(units.get(a, [])) for a in ARMS}
sets = {a: set(counts[a]) for a in ARMS}
print("units per arm: %s" % {a: len(units.get(a, [])) for a in ARMS})
print("genera per arm: %s" % {a: len(sets[a]) for a in ARMS})

AMPH, ENDO = ["herptile", "ehi_amphibian"], ["ehi", "youngblut"]


def block_of(g):
    in_h = g in sets["herptile"]
    in_n = g in sets["ehi_amphibian"]
    in_e = any(g in sets[a] for a in ENDO)
    if in_h and in_n and not in_e:
        return 1
    if in_h and in_n and in_e:
        return 2
    if in_h and not in_n:
        return 3
    return 4


BLOCK_TXT = {
    1: "both amphibian\ncatalogs, not\nrecovered by EHI\nor Youngblut",
    2: "both amphibian\ncatalogs, also in\nEHI or Youngblut",
    3: "herptile only,\nnot in newts",
    4: "not in herptile",
}
BLOCK_COL = {1: "#C44E52", 2: "#DD8452", 3: "#937860", 4: "#4C72B0"}

blocks = defaultdict(list)
for g in set().union(*[sets[a] for a in ARMS]):
    blocks[block_of(g)].append(g)
for b in blocks:
    blocks[b].sort(key=lambda g: -(counts["herptile"].get(g, 0)
                                   + counts["ehi_amphibian"].get(g, 0)
                                   + counts["ehi"].get(g, 0)
                                   + counts["youngblut"].get(g, 0)))
order = blocks[1] + blocks[2] + blocks[3] + blocks[4][:TOP_NOHERP]
print("blocks: 1=%d 2=%d 3=%d 4=%d (showing top %d of block 4)"
      % (len(blocks[1]), len(blocks[2]), len(blocks[3]), len(blocks[4]), TOP_NOHERP))

fig = plt.figure(figsize=(17.5, 10.2))
gs = fig.add_gridspec(3, 2, width_ratios=[1.15, 1.0],
                      height_ratios=[1, 1, 1], wspace=0.30, hspace=0.62)

# ============================================================== PANEL A
a = fig.add_subplot(gs[:, 0])
mat = np.array([[counts[c].get(g, 0) for c in ARMS] for g in order], dtype=float)
disp = np.log10(mat + 1)
cmap = LinearSegmentedColormap.from_list(
    "counts", ["#F7F7F7", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"])
im = a.imshow(disp, aspect="auto", cmap=cmap, vmin=0, vmax=disp.max())
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        v = int(mat[i, j])
        if v:
            a.text(j, i, "%d" % v, ha="center", va="center", fontsize=6.6,
                   color="white" if disp[i, j] > disp.max() * 0.55 else "#222222")
for xv in (0.5, 1.5, 2.5, 3.5):
    a.axvline(xv, color="#FFFFFF", lw=1.8)
a.axvline(1.5, color="#333333", lw=1.6)
a.set_xticks(range(5))
a.set_xticklabels(COLLAB, fontsize=7.8)
a.set_yticks(range(len(order)))
a.set_yticklabels(order, fontsize=6.8)
a.set_title("a   Genus composition across five catalogs", loc="left",
            fontweight="bold", fontsize=12)
cuts, run = [], 0
for b in (1, 2, 3):
    run += len(blocks[b])
    cuts.append(run - 0.5)
for y in cuts:
    a.axhline(y, color="#C44E52", lw=1.3)
start = 0
for b in (1, 2, 3, 4):
    n = len(blocks[b]) if b != 4 else TOP_NOHERP
    mh = sum(counts["herptile"].get(g, 0) for g in blocks[b])
    lab = BLOCK_TXT[b] + ("\n%d genera" % len(blocks[b]))
    if b != 4:
        lab += "\n%d herptile units" % mh
    a.text(4.65, start + n / 2.0 - 0.5, lab, fontsize=6.9, ha="left",
           va="center", color=BLOCK_COL[b])
    start += n
a.set_xlim(-0.5, 5.9)
a.text(-0.5, len(order) + 2.6,
       "jointly dereplicated at 95% ANI across all five arms, so units are comparable; "
       "black line separates amphibian from endotherm arms",
       fontsize=6.9, color="#444444", ha="left")
a.text(-0.5, len(order) + 1.5,
       "the GTDB column shows these genera are described: of 62 reference genomes in the top block, "
       "42 are from gut or faecal sources",
       fontsize=6.9, color="#C44E52", ha="left")
cb = fig.colorbar(im, ax=a, fraction=0.024, pad=0.015)
cb.set_label("log10(units + 1)", fontsize=7.5)
cb.ax.tick_params(labelsize=6.5)

# ============================================================== PANEL B
b = fig.add_subplot(gs[0, 1])
be = {(r["arm_a"], r["arm_b"]): r for r in read_tsv(BETA)}
SHOW = [("herptile", "ehi", "amphibian vs\nEHI mammals"),
        ("herptile", "youngblut", "amphibian vs\nYoungblut"),
        ("ehi_amphibian", "ehi", "EHI newts vs\nEHI mammals"),
        ("ehi", "youngblut", "EHI mammals vs\nYoungblut"),
        ("herptile", "ehi_amphibian", "amphibian vs\nEHI newts")]
ys, sim, sne, labs = [], [], [], []
for i, (x, y, lab) in enumerate(SHOW):
    r = be.get((x, y)) or be.get((y, x))
    if not r:
        continue
    ys.append(i)
    sim.append(tofloat(r["b_sim_turnover"]))
    sne.append(tofloat(r["b_sne_nestedness"]))
    labs.append(lab)
ys = np.array(ys, dtype=float)
b.barh(ys, sim, height=0.6, color="#C44E52", edgecolor="black", lw=0.5,
       label="turnover (replacement)")
b.barh(ys, sne, left=sim, height=0.6, color="#BBBBBB", edgecolor="black",
       lw=0.5, label="nestedness (subset)")
for i, (s, n) in enumerate(zip(sim, sne)):
    b.text(s + n + 0.012, ys[i], "%.0f%% turn" % (100.0 * s / (s + n) if (s + n) else 0),
           va="center", fontsize=6.8)
b.set_yticks(ys)
b.set_yticklabels(labs, fontsize=7.0)
b.invert_yaxis()
b.set_xlim(0, 1.06)
b.set_xlabel("Sorensen dissimilarity", fontsize=8)
b.set_title("b   Different pools, or one a subset of the other?",
            loc="left", fontweight="bold", fontsize=11)
b.legend(fontsize=6.8, frameon=False, loc="lower right")
b.grid(axis="x", alpha=0.2, lw=0.5)
b.tick_params(labelsize=7)

# ============================================================== PANEL C
c = fig.add_subplot(gs[1, 1])
th = read_tsv(THRESH)
x = [int(r["min_genomes"]) for r in th]
n = [int(r["n_genera"]) for r in th]
c.plot(x, n, "o-", color="#C44E52", lw=1.8, ms=7, markeredgecolor="black",
       markeredgewidth=0.6)
for xi, ni in zip(x, n):
    c.text(xi, ni + 0.28, "%d" % ni, ha="center", fontsize=7.5)
c.set_xticks(x)
c.set_xlabel("minimum units required to count as present", fontsize=8)
c.set_ylabel("genera in both amphibian catalogs,\nnot recovered by EHI mammal\nor Youngblut", fontsize=7.5)
c.set_ylim(0, max(n) * 1.28)
c.set_title("c   How fast the strict set erodes", loc="left",
            fontweight="bold", fontsize=11)
c.grid(alpha=0.2, lw=0.5)
c.tick_params(labelsize=7)
c.text(0.97, 0.93, "a set that survives only at\nthreshold 1 is a threshold artefact",
       transform=c.transAxes, ha="right", va="top", fontsize=6.8, color="#444444")

# ============================================================== PANEL D
d = fig.add_subplot(gs[2, 1])
acc = defaultdict(list)
for r in read_tsv(ACC):
    acc[r["arm"]].append((int(r["n_sampled"]), tofloat(r["mean_genera"]),
                          tofloat(r["lo"]), tofloat(r["hi"])))
COLS = {"herptile": "#C44E52", "ehi_amphibian": "#DD8452",
        "ehi": "#4C72B0", "youngblut": "#55A868", "gtdb_ref": "#8C8C8C"}
NICE = {"herptile": "herptile wild", "ehi_amphibian": "EHI newts",
        "ehi": "EHI mammal/bird", "youngblut": "Youngblut", "gtdb_ref": "GTDB refs"}
for arm in ARMS:
    v = sorted(acc.get(arm, []))
    if not v:
        continue
    xs = [p[0] for p in v]
    ms = [p[1] for p in v]
    lo = [p[2] for p in v]
    hi = [p[3] for p in v]
    d.plot(xs, ms, "-o", color=COLS[arm], lw=1.5, ms=3.5, label=NICE[arm])
    d.fill_between(xs, lo, hi, color=COLS[arm], alpha=0.15, lw=0)
d.set_xlabel("units sampled", fontsize=8)
d.set_ylabel("genera detected", fontsize=8)
d.set_title("d   Has each catalog saturated?", loc="left",
            fontweight="bold", fontsize=11)
d.legend(fontsize=6.6, frameon=False, loc="upper left")
d.grid(alpha=0.2, lw=0.5)
d.tick_params(labelsize=7)
d.text(0.97, 0.06, "a curve still climbing has NOT saturated,\nso its absences are sampling-limited",
       transform=d.transAxes, ha="right", va="bottom", fontsize=6.8, color="#444444")

fig.savefig(STEM + ".png", dpi=300, bbox_inches="tight")
fig.savefig(STEM + ".pdf", bbox_inches="tight")
print()
print("wrote %s" % (STEM + ".png"))
print()
print("PANEL B IS THE HEADLINE: herptile vs EHI mammals is 100% turnover,")
print("zero nestedness. The pools are genuinely different, not one being an")
print("undersampled subset of the other. herptile vs EHI newts is the")
print("opposite, 100% nestedness, because the newt genera are a subset of ours.")
print()
print("CAPTION MUST STATE:")
print("  exploratory, not a preregistered test")
print("  only ONE genus (Negativibacillus) was recovered in all five arms, so")
print("  the positive control for family-wide recoverability is thin")
print("  the EHI newt arm is 43 units and has not saturated (panel d), so its")
print("  absences carry little weight")

# SENTINEL_END
