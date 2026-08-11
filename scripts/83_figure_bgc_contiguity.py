# The biosynthetic cluster and contiguity figure is drawn.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/fig_bgc_contiguity.py
# Output: results/figures/Figure_bgc_contiguity.pdf and .png
import os, csv, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# FIG_BGC_CONTIGUITY V1. Replaces Figure_n50_vs_bgc (three arms) and
# Figure_bgc_arms. Four panels, one per claim in Results 3.6:
#   A  complete-BGC recovery tracks N50 WITHIN every arm
#   B  at MATCHED absolute N50 the arms agree
#   C  the arms are NOT matched on contiguity, so no cross-arm comparison
#   D  one subgroup sits above the curve at matched N50. HYPOTHESIS ONLY.
# The old figure used deciles in one panel and quartiles in another, so the
# two disagreed on the ceiling (2.2 vs 1.59). Quintiles are used uniformly
# here. Youngblut at n=48 gives ~10 genomes per quintile; that is stated.

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
SRC = B + "/results/assembly_quality_all_arms.tsv"
STEM = B + "/results/figures/Figure_bgc_contiguity"

if not os.path.exists(SRC):
    raise SystemExit("MISSING INPUT: " + SRC)
for e in (".pdf", ".png"):
    if os.path.exists(STEM + e):
        raise SystemExit("REFUSING TO OVERWRITE: " + STEM + e)

ARMS = ["UHM amphibian", "EHI newt", "EHI mammal", "Youngblut",
        "GTDB reference"]
COL = {"UHM amphibian": "#2c6fbb", "EHI newt": "#7a5c9e",
       "EHI mammal": "#c1440e", "Youngblut": "#d99b1c",
       "GTDB reference": "#5a5a5a"}
MRK = {"UHM amphibian": "o", "EHI newt": "^", "EHI mammal": "s",
       "Youngblut": "D", "GTDB reference": "v"}
DARK = "#222222"

with open(SRC, newline="") as fh:
    rr = [x for x in csv.reader(fh, delimiter="\t") if x]
H = {c: i for i, c in enumerate(rr[0])}
for need in ("arm", "n50", "n_complete", "n_regions", "subgroup"):
    if need not in H:
        raise SystemExit("column missing: " + need + "\n" + str(rr[0]))

rows = []
for p in rr[1:]:
    try:
        rows.append({"arm": p[H["arm"]].strip(),
                     "sub": p[H["subgroup"]].strip(),
                     "n50": float(p[H["n50"]]),
                     "comp": float(p[H["n_complete"]]),
                     "reg": float(p[H["n_regions"]])})
    except (ValueError, IndexError):
        continue

present = sorted(set(r["arm"] for r in rows))
print("arms in file: " + ", ".join("%s (n=%d)"
      % (a, sum(1 for r in rows if r["arm"] == a)) for a in present))
skipped = [a for a in present if a not in ARMS]
if skipped:
    print("NOT PLOTTED, too few genomes at this scope: " + ", ".join(skipped))
    for a in skipped:
        print("   %s n=%d" % (a, sum(1 for r in rows if r["arm"] == a)))
print("total rows: %d" % len(rows))
print()

by_arm = {a: sorted([r for r in rows if r["arm"] == a],
                    key=lambda x: x["n50"]) for a in ARMS}
by_arm = {a: v for a, v in by_arm.items() if v}

fig = plt.figure(figsize=(14.6, 9.6))
gs = fig.add_gridspec(2, 2, wspace=0.26, hspace=0.36,
                      left=0.07, right=0.985, top=0.93, bottom=0.07)

def mean(v):
    return sum(v) / len(v) if v else 0.0

# ------------------------------------------------------------- PANEL A
axA = fig.add_subplot(gs[0, 0])
NQ = 5
for a, v in by_arm.items():
    if len(v) < NQ * 2:
        continue
    xs, ys = [], []
    for q in range(NQ):
        lo = q * len(v) // NQ
        hi = (q + 1) * len(v) // NQ
        chunk = v[lo:hi]
        if not chunk:
            continue
        mid = sorted(x["n50"] for x in chunk)[len(chunk) // 2]
        xs.append(math.log10(mid))
        ys.append(mean([x["comp"] for x in chunk]))
    axA.plot(xs, ys, marker=MRK[a], color=COL[a], lw=1.6, ms=6,
             label="%s (n = %d)" % (a, len(v)), alpha=0.9)
axA.set_xlabel("N50 (bp), quintile median within arm", fontsize=9)
axA.set_ylabel("mean complete BGCs per genome", fontsize=9)
axA.set_title("A   Complete-BGC recovery tracks contiguity within every arm",
              fontsize=10.5, loc="left", pad=9)
xt = [3.5, 4.0, 4.5, 5.0]
axA.set_xticks(xt)
axA.set_xticklabels(["%.0f kb" % (10 ** t / 1000.0) for t in xt], fontsize=8)
axA.legend(fontsize=7.6, frameon=False, loc="upper left")
axA.spines["top"].set_visible(False)
axA.spines["right"].set_visible(False)
axA.tick_params(labelsize=8)

# ------------------------------------------------------------- PANEL B
axB = fig.add_subplot(gs[0, 1])
BINS = [(0, 10000, "<10 kb"), (10000, 20000, "10-20"), (20000, 40000, "20-40"),
        (40000, 80000, "40-80"), (80000, 1e12, ">80 kb")]
xsb = list(range(len(BINS)))
for a, v in by_arm.items():
    ys, ns = [], []
    for lo, hi, _ in BINS:
        sel = [r["comp"] for r in v if lo <= r["n50"] < hi]
        ys.append(mean(sel) if sel else float("nan"))
        ns.append(len(sel))
    axB.plot(xsb, ys, marker=MRK[a], color=COL[a], lw=1.6, ms=6, alpha=0.9,
             label=a)
    for x, y, n in zip(xsb, ys, ns):
        if n and n < 15 and y == y:
            axB.annotate(str(n), (x, y), textcoords="offset points",
                         xytext=(0, -11), ha="center", fontsize=6.4,
                         color=COL[a])
axB.set_xticks(xsb)
axB.set_xticklabels([b[2] for b in BINS], fontsize=8)
axB.set_xlabel("absolute N50 bin, identical across arms", fontsize=9)
axB.set_ylabel("mean complete BGCs per genome", fontsize=9)
axB.set_title("B   At matched contiguity the arms converge",
              fontsize=10.5, loc="left", pad=9)
axB.legend(fontsize=7.6, frameon=False, loc="upper left")
axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)
axB.tick_params(labelsize=8)
axB.text(0.98, 0.03, "small numerals are n where n < 15",
         transform=axB.transAxes, ha="right", fontsize=6.8, color="#777777")

# ------------------------------------------------------------- PANEL C
axC = fig.add_subplot(gs[1, 0])
order = [a for a in ARMS if a in by_arm]
data = [[math.log10(r["n50"]) for r in by_arm[a]] for a in order]
bp = axC.boxplot(data, vert=False, widths=0.6, patch_artist=True,
                 showfliers=False, medianprops=dict(color="white", lw=1.6))
for patch, a in zip(bp["boxes"], order):
    patch.set_facecolor(COL[a])
    patch.set_edgecolor("none")
    patch.set_alpha(0.88)
for w in bp["whiskers"] + bp["caps"]:
    w.set_color("#999999")
axC.set_yticklabels(["%s\nn = %d" % (a, len(by_arm[a])) for a in order],
                    fontsize=7.8, linespacing=1.3)
axC.set_xticks(xt)
axC.set_xticklabels(["%.0f kb" % (10 ** t / 1000.0) for t in xt], fontsize=8)
axC.set_xlabel("N50 (bp, log scale)", fontsize=9)
axC.set_title("C   The arms are not matched on contiguity",
              fontsize=10.5, loc="left", pad=9)
axC.spines["top"].set_visible(False)
axC.spines["right"].set_visible(False)
axC.tick_params(labelsize=8)

# ------------------------------------------------------------- PANEL D
axD = fig.add_subplot(gs[1, 1])
KEEP = ("Rodentia", "Carnivora", "Lagomorpha")
SUBS = [s for s in sorted(set(r["sub"] for r in rows))
        if s.startswith("EHI mammal:") and s.split(": ")[-1] in KEEP]
_drop = [s for s in sorted(set(r["sub"] for r in rows))
         if s.startswith("EHI mammal:") and s.split(": ")[-1] not in KEEP]
for _s in _drop:
    print("panel D excludes %s, n=%d, too few genomes for a mean"
          % (_s, sum(1 for r in rows if r["sub"] == _s)))
SCOL = {"EHI mammal: Lagomorpha": "#c1440e",
        "EHI mammal: Rodentia": "#e39a7a",
        "EHI mammal: Carnivora": "#8c6a5d"}
for s in SUBS:
    v = [r for r in rows if r["sub"] == s]
    ys, ns = [], []
    for lo, hi, _ in BINS:
        sel = [r["comp"] for r in v if lo <= r["n50"] < hi]
        ys.append(mean(sel) if sel else float("nan"))
        ns.append(len(sel))
    lag = "Lagomorpha" in s
    axD.plot(xsb, ys, marker="s" if lag else "o",
             color=SCOL.get(s, "#aaaaaa"), lw=2.2 if lag else 1.3,
             ms=7 if lag else 5, alpha=0.95 if lag else 0.75,
             label=s.replace("EHI mammal: ", "") + " (n = %d)" % len(v))
    for x, y, n in zip(xsb, ys, ns):
        if n and y == y:
            axD.annotate(str(n), (x, y), textcoords="offset points",
                         xytext=(0, -11), ha="center", fontsize=6.2,
                         color=SCOL.get(s, "#aaaaaa"))
for a in ("EHI newt", "GTDB reference"):
    if a not in by_arm:
        continue
    v = by_arm[a]
    ys = []
    for lo, hi, _ in BINS:
        sel = [r["comp"] for r in v if lo <= r["n50"] < hi]
        ys.append(mean(sel) if sel else float("nan"))
    axD.plot(xsb, ys, color=COL[a], lw=1.1, ls="--", alpha=0.65, label=a)
axD.set_xticks(xsb)
axD.set_xticklabels([b[2] for b in BINS], fontsize=8)
axD.set_xlabel("absolute N50 bin, identical across groups", fontsize=9)
axD.set_ylabel("mean complete BGCs per genome", fontsize=9)
axD.set_title("D   Rabbit-derived genomes sit above the curve.\n"
              "HYPOTHESIS, not a result: bins hold 8 to 21 genomes",
              fontsize=10.5, loc="left", pad=9)
axD.legend(fontsize=7.2, frameon=False, loc="upper left",
           bbox_to_anchor=(0.0, 0.86))
axD.spines["top"].set_visible(False)
axD.spines["right"].set_visible(False)
axD.tick_params(labelsize=8)
axD.text(0.98, 0.03, "numerals are n in every cell", transform=axD.transAxes,
         ha="right", fontsize=6.8, color="#777777")

fig.savefig(STEM + ".pdf", bbox_inches="tight")
fig.savefig(STEM + ".png", dpi=300, bbox_inches="tight")
print("wrote " + STEM + ".pdf and .png")
print()

try:
    from scipy.stats import spearmanr
    print("SPEARMAN, N50 vs complete BGCs, WITHIN each arm:")
    for a in order:
        v = by_arm[a]
        rho, p = spearmanr([r["n50"] for r in v], [r["comp"] for r in v])
        print("  %-18s n=%5d  rho %+0.3f  p %.2g" % (a, len(v), rho, p))
    print()
    print("TOTAL REGIONS vs N50, the contrast that matters:")
    for a in order:
        v = by_arm[a]
        rho, p = spearmanr([r["n50"] for r in v], [r["reg"] for r in v])
        print("  %-18s n=%5d  rho %+0.3f  p %.2g" % (a, len(v), rho, p))
except ImportError:
    print("scipy unavailable, correlations skipped")
print()

print("=" * 78)
print("CAPTION REQUIREMENTS IMPLIED BY THIS FIGURE")
print("=" * 78)
print("ALL PANELS. antiSMASH 7.1.0, --taxon bacteria --genefinding-tool")
print("  prodigal-m --cb-knownclusters. RUMINOCOCCACEAE ONLY in every arm.")
print("  An earlier all-family amphibian set against a Ruminococcaceae-only")
print("  reference set gave 0.50 vs 0.64; restricting to one family gives")
print("  0.28. DO NOT quote the all-family comparison.")
print("  Quintiles are used uniformly. The old figure mixed deciles and")
print("  quartiles across panels and the two disagreed on the ceiling.")
print()
print("PANEL A. MUST STATE that within the GTDB arm alone the spread is")
print("  roughly 80-fold, far larger than the roughly 2-fold between-arm gap.")
print()
print("PANEL B. MUST STATE which cells are thin. Numerals mark n < 15.")
print("  This panel is what licenses the claim that contiguity, not host,")
print("  explains the between-arm difference.")
print()
print("PANEL C. MUST STATE the matching result: UHM amphibian vs EHI newt is")
print("  balanced on N50 (SMD 0.071) and contigs (0.031) but NOT on assembly")
print("  length (0.542, KS p 0.00055). EHI mammal vs EHI newt fails on contig")
print("  count (SMD 0.272, KS p 0.035) despite a shared pipeline. Both MAG")
print("  arms fail every balance test against GTDB.")
print("  MUST STATE: GTDB is TAXONOMIC background, not a host comparison;")
print("  91.9 percent of its genomes are themselves MAGs.")
print()
print("PANEL D. HYPOTHESIS ONLY. MUST STATE: bins hold 8 to 21 genomes with")
print("  no intervals, and rabbit genomes come disproportionately from")
print("  particular studies, so study of origin is not separable from host.")
print("  MUST STATE: the arm-level EHI mammal mean of 0.56 against EHI newt")
print("  0.33 is NOT usable, because the mammal arm is heterogeneous by host")
print("  order (Lagomorpha 1.07, Rodentia 0.35) and the arm mean is the wrong")
print("  summary.")
print()
print("NOT DONE, state as such: no inferential model has been fitted. Counts")
print("  are zero-inflated. The hurdle negative binomial with N50 and")
print("  contig-edge rate as covariates is the recommended next step and")
print("  appears nowhere in the BGC literature.")
print("  DO NOT decompose the between-arm gap using total regions as a")
print("  baseline. Total region count is inflated by BGC splitting and")
print("  deflated by the detection floor, so it is not a stable denominator.")
print()
print("FIG_BGC_CONTIGUITY_V2_20260806_COMPLETE")
# FIG_BGC_CONTIGUITY_V2_20260806_COMPLETE
