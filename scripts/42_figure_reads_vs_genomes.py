# The reads versus genomes figure is drawn.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/fig_reads_vs_genomes.py
# Output: results/figures/Figure_reads_vs_genomes.pdf and .png
import os, csv, math, random
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# FIG_READS_VS_GENOMES V2. V1 fixes: label offsets were hand-set and four
# pointed at the wrong marker (Fimivivens, Limiplasma, Anaerostipes, and
# Ruthenibacterium overlapping Intestinimonas); the "not detected" note
# collided with the legend; the y-axis label implied y=0 meant absent from
# the catalog when it means captive-derived genomes only; and literal %%
# leaked into the printed caption text.

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
PHY = B + "/results/bracken_wild_phylum.tsv"
REC = B + "/results/bracken_reciprocal_recovery_v2.tsv"
MAP = B + "/results/gtdb_ncbi_reciprocal.tsv"
SGB = B + "/data/sgb_manifest.tsv"
STEM = B + "/results/figures/Figure_reads_vs_genomes"

for p in (PHY, REC, MAP, SGB):
    if not os.path.exists(p):
        raise SystemExit("MISSING INPUT: " + p)
for e in (".pdf", ".png"):
    if os.path.exists(STEM + e):
        raise SystemExit("REFUSING TO OVERWRITE: " + STEM + e)

def read_tsv(path):
    with open(path, newline="") as fh:
        r = [x for x in csv.reader(fh, delimiter="\t") if x]
    return r[0], r[1:]

WILD = "#2c6fbb"
CAPT = "#c1440e"
NORD = "#7a5c9e"
GREY = "#9a9a9a"
DARK = "#222222"

fig = plt.figure(figsize=(15.6, 5.6))
gs = fig.add_gridspec(1, 3, width_ratios=[0.92, 0.85, 1.75], wspace=0.40,
                      left=0.06, right=0.985, top=0.84, bottom=0.16)

# ------------------------------------------------------------- PANEL A
axA = fig.add_subplot(gs[0, 0])
h, rows = read_tsv(PHY)
I = {c: i for i, c in enumerate(h)}
ph = []
for p in rows:
    try:
        ph.append((p[I["phylum"]], float(p[I["mean_fraction_wild"]])))
    except (ValueError, KeyError, IndexError):
        continue
ph.sort(key=lambda x: -x[1])
ph = ph[:8]
ys = list(range(len(ph)))[::-1]
for y, (name, v) in zip(ys, ph):
    c = WILD if name == "Bacillota" else GREY
    axA.barh(y, v, height=0.68, color=c, edgecolor="none")
    axA.text(v + 0.008, y, "%.3f" % v, va="center", fontsize=8,
             color=DARK if name == "Bacillota" else "#666666",
             fontweight="bold" if name == "Bacillota" else "normal")
axA.set_yticks(ys)
axA.set_yticklabels([n for n, _ in ph], fontsize=8.4)
axA.set_xlabel("mean read fraction, 44 wild samples", fontsize=9)
axA.set_xlim(0, max(v for _, v in ph) * 1.28)
axA.set_title("A   Bacillota is under a quarter of\nthese communities",
              fontsize=10, loc="left", pad=9)
axA.spines["top"].set_visible(False)
axA.spines["right"].set_visible(False)
axA.tick_params(axis="both", labelsize=8)

# ------------------------------------------------------------- PANEL B
axB = fig.add_subplot(gs[0, 1])
h, rows = read_tsv(SGB)
J = {c: i for i, c in enumerate(h)}
wild_g = set()
for p in rows:
    g = p[J["genus"]].strip()
    if g and g.upper() != "UNASSIGNED" and p[J["has_wild"]].strip() == "yes":
        wild_g.add(g)

h, rows = read_tsv(MAP)
K = {c: i for i, c in enumerate(h)}
verdict = {}
for p in rows:
    verdict[p[K["gtdb_genus"]].strip()] = p[K["verdict_t70"]].strip()

CATS = [("reciprocal", "testable", WILD),
        ("ncbi_name_is_shared", "NCBI name shared\nacross GTDB genera", GREY),
        ("no_ncbi_counterpart", "no NCBI name on any\nreference genome", GREY),
        ("both_directions_fail", "both directions fail", GREY),
        ("forward_ambiguous", "forward ambiguous", GREY)]
cnt = {k: 0 for k, _, _ in CATS}
for g in wild_g:
    v = verdict.get(g, "no_ncbi_counterpart")
    if v in ("not_in_gtdb", "not_in_reciprocal_file"):
        v = "no_ncbi_counterpart"
    if v in cnt:
        cnt[v] += 1
TOT = len(wild_g)

ys = list(range(len(CATS)))[::-1]
for y, (k, lab, c) in zip(ys, CATS):
    axB.barh(y, cnt[k], height=0.66, color=c, edgecolor="none")
    axB.text(cnt[k] + 0.8, y, str(cnt[k]), va="center", fontsize=8.4,
             color=DARK if k == "reciprocal" else "#666666",
             fontweight="bold" if k == "reciprocal" else "normal")
axB.set_yticks(ys)
axB.set_yticklabels([l for _, l, _ in CATS], fontsize=7.6, linespacing=1.3)
axB.set_xlim(0, max(cnt.values()) * 1.25)
axB.set_xlabel("wild catalog genera (n = %d)" % TOT, fontsize=9)
axB.set_title("B   Only %d of %d wild genera can be\ncompared with read profiles"
              % (cnt["reciprocal"], TOT), fontsize=10, loc="left", pad=9)
axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)
axB.tick_params(axis="both", labelsize=8)

# ------------------------------------------------------------- PANEL C
axC = fig.add_subplot(gs[0, 2])
h, rows = read_tsv(REC)
L = {c: i for i, c in enumerate(h)}
recs = []
for p in rows:
    try:
        recs.append({"g": p[L["gtdb_genus"]].strip(),
                     "b": p[L["bucket"]].strip(),
                     "f": float(p[L["mean_fraction_wild"]]),
                     "w": int(p[L["wild_sgbs"]])})
    except (ValueError, KeyError, IndexError):
        continue

KEEP = ("reads_and_wild_genome", "wild_genome_no_reads")
recs = [r for r in recs if r["b"] in KEEP]
print("panel C restricted to %d wild catalog genera" % len(recs))
nz = [r["f"] for r in recs if r["f"] > 0]
LO, HI = math.log10(min(nz)), math.log10(max(nz))
XZERO = LO - 0.62

def X(f):
    return math.log10(f) if f > 0 else XZERO

def Y(n):
    return math.log10(1 + n)

YTOP = Y(42) * 1.16
axC.axvspan(XZERO - 0.30, XZERO + 0.22, color="#f2f2f2", zorder=0)
axC.axhline(Y(0), color="#e0e0e0", lw=1.0, zorder=0)

STYLE = [("reads_and_wild_genome", WILD, "o", 62, "reads and a wild genome"),
         ("wild_genome_no_reads", NORD, "^", 60, "wild genome, no reads")]

random.seed(7)
placed = []
for key, c, m, s, lab in STYLE:
    sel = [r for r in recs if r["b"] == key]
    if not sel:
        continue
    xs, ys2 = [], []
    for r in sel:
        jx = random.uniform(-0.11, 0.11) if r["f"] == 0 else 0.0
        jy = random.uniform(-0.022, 0.022) if r["w"] == 0 else 0.0
        px, py = X(r["f"]) + jx, Y(r["w"]) + jy
        xs.append(px)
        ys2.append(py)
        placed.append((r["g"], px, py))
    axC.scatter(xs, ys2, s=s, marker=m, color=c, alpha=0.85, zorder=3,
                edgecolor="white", linewidth=0.7,
                label=lab + "  (n = %d)" % len(sel))

# Labels placed against the ACTUAL plotted coordinate, with collision
# avoidance, instead of hand-set offsets that drifted onto wrong markers.
NAMED = ["Anaerotruncus", "Intestinimonas", "Enterocloster", "Oscillibacter",
         "Ruthenibacterium", "Faecalibacterium", "Anaerostipes",
         "Angelakisella", "Evtepia", "Fimivivens"]
coord = {g: (x, y) for g, x, y in placed}
XSPAN = (HI + 0.30) - (XZERO - 0.36)
taken = []

def collides(x, y):
    for tx, ty in taken:
        near_zero = abs(y) < 0.09 * YTOP and abs(ty) < 0.09 * YTOP
        ytol = 0.150 * YTOP if near_zero else 0.075 * YTOP
        if abs(x - tx) < 0.34 * XSPAN / 6.0 and abs(y - ty) < ytol:
            return True
    return False

for g in NAMED:
    if g not in coord:
        continue
    x, y = coord[g]
    right = x < (HI - 0.35)
    for dy in (0.045, -0.055, 0.105, -0.115, 0.165, -0.175):
        lx = x + (0.085 if right else -0.085)
        ly = y + dy * YTOP
        if not collides(lx, ly) and 0 < ly < YTOP:
            taken.append((lx, ly))
            axC.annotate(g, (x, y), xytext=(lx, ly), fontsize=7.8, color=DARK,
                         ha="left" if right else "right", va="center",
                         arrowprops=dict(arrowstyle="-", lw=0.55,
                                         color="#aaaaaa", shrinkA=0,
                                         shrinkB=3))
            break

ticks, labs = [XZERO], ["none"]
t = math.ceil(LO * 2) / 2.0
while t <= HI + 0.01:
    ticks.append(t)
    labs.append("%.4f" % (10 ** t))
    t += 0.5
axC.set_xticks(ticks)
axC.set_xticklabels(labs, fontsize=7.8)
axC.set_xlim(XZERO - 0.36, HI + 0.30)

yt = [0, 1, 2, 5, 10, 20, 42]
axC.set_yticks([Y(v) for v in yt])
axC.set_yticklabels([str(v) for v in yt], fontsize=8)
axC.set_ylim(-0.30, YTOP)

axC.set_xlabel("mean read fraction across 44 wild samples", fontsize=9)
axC.set_ylabel("wild SGBs recovered", 
               fontsize=9, linespacing=1.4)
axC.set_title("C   Read abundance is a poor predictor of genome recovery", fontsize=10, loc="left", pad=9)
axC.spines["top"].set_visible(False)
axC.spines["right"].set_visible(False)
axC.legend(loc="upper left", fontsize=7.6, frameon=False,
           handletextpad=0.5, borderpad=0.2)
axC.text(XZERO, -0.275, "no reads assigned", ha="center", va="bottom",
         fontsize=7.2, color="#777777")

fig.savefig(STEM + ".pdf", bbox_inches="tight")
fig.savefig(STEM + ".png", dpi=300, bbox_inches="tight")
print("wrote " + STEM + ".pdf and .png")

CAP = """
==============================================================================
CAPTION REQUIREMENTS IMPLIED BY THIS FIGURE
==============================================================================
PANEL A. Bracken pluspf_20251015, phylum level, 44 wild samples.
  MUST STATE: NCBI Bacillota is BROADER than GTDB Bacillota_A, so 0.223 is
  an UPPER BOUND on the clade this catalog covers.
  MUST STATE: reads were CLASSIFIED to a taxon. That is not the same as
  demonstrating the taxon was present.

PANEL B. Testability of the {tot} wild catalog genera.
  MUST STATE the criterion: a GTDB-NCBI genus pair counts as comparable only
  if the GTDB genus is at least 70 percent labelled with that NCBI genus AND
  genomes carrying that NCBI name are at least 70 percent inside that GTDB
  genus. A forward-only rule is NOT sufficient: GTDB Merdicola maps onto NCBI
  Clostridium, which spans 213 GTDB genera.
  MUST STATE: the 70 percent threshold was fixed before the join was run.
  MUST STATE: "no NCBI name on any reference genome" does NOT mean the genus
  is absent from GTDB.
  MUST STATE: do NOT collapse the untestable categories into "not detected".
  {untest} of {tot} genera were never testable.

PANEL C. The {rec} testable wild catalog genera, two outcomes.
{buckets}
  Points in the shaded column had NO reads assigned in ANY of the 44 wild
  samples. Points at y = 0 were not recovered as wild genomes.
  Both are real categories, not missing data.
  Jitter is applied at zero positions only.
  MUST STATE the key contrast: Faecalibacterium reads were about TWICE as
  abundant as Anaerotruncus reads in the same samples, yet 42 wild
  Anaerotruncus genomes were recovered and no wild Faecalibacterium genome.
  MUST STATE the reverse: Angelakisella has 21 wild SGBs and zero read
  assignments, expected where the classification database holds no
  sufficiently close reference genome.

WHOLE FIGURE. Catalog side is ALL 1,171 SGBs at genus level with NO family
  filter. A family-level comparison is impossible here: 96.2 percent of GTDB
  Ruminococcaceae reference genomes carry NCBI family Oscillospiraceae, as do
  94.7 percent of GTDB Oscillospiraceae.

DO NOT report the 8.6 percent "share of Bacillota reads" figure. Its
  numerator covers only the 14 recovered genera while its denominator is all
  Bacillota reads, so it is not a recovery fraction.
"""
bl = "\n".join("    %-38s %d" % (lab, sum(1 for r in recs if r["b"] == k))
               for k, c, m, s, lab in STYLE)
print(CAP.format(tot=TOT, untest=TOT - cnt["reciprocal"],
                 rec=cnt["reciprocal"], buckets=bl))
print("FIG_READS_VS_GENOMES_V3_20260806_COMPLETE")
# FIG_READS_VS_GENOMES_V3_20260806_COMPLETE
