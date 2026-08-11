# The gene content figure is drawn.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/fig_gene_content.py
# Output: results/figures/Figure_gene_content.pdf and .png
import os, csv, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# FIG_GENE_CONTENT V2. Three panels for Results 3.5.
# V1 errors fixed: panel A title was a hard-coded literal that disagreed with
# the data; WHOLE_TREE was plotted as if it were a genus; the highlight rule
# keyed on the verdict string, which calls underpowered non-significant tests
# "testable". The gate is now the independent transition count.

DROP_ENDO = False
MIN_TRANSITIONS = 5

BASE = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
INTER = os.path.join(BASE, "results", "genus_interleaving.tsv")
HAPPI = os.path.join(BASE, "work", "focal_genus_pangenome", "matrices",
                     "happi_results_og_bacteria.tsv")
HMM = os.path.join(BASE, "work", "focal_genus_pangenome", "pyrg_control", "hmm",
                   "hmm_genome_model.tsv")
Q = os.path.join(BASE, "work", "focal_genus_pangenome", "checkm2_out",
                 "quality_report.tsv")
OUTD = os.path.join(BASE, "results", "figures")
STEM = os.path.join(OUTD, "Figure_gene_content")

for p in (INTER, HAPPI, HMM):
    if not os.path.exists(p):
        raise SystemExit("MISSING INPUT: " + p)
for ext in (".pdf", ".png"):
    if os.path.exists(STEM + ext):
        raise SystemExit("REFUSING TO OVERWRITE: " + STEM + ext)

def read_tsv(path):
    with open(path, newline="") as fh:
        r = [x for x in csv.reader(fh, delimiter="\t") if x]
    return r[0], r[1:]

AMPH = "#2c6fbb"
REF = "#c1440e"
GREY = "#9a9a9a"
DARK = "#222222"

fig = plt.figure(figsize=(15.4, 5.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.2, 1.4], wspace=0.46,
                      left=0.055, right=0.985, top=0.85, bottom=0.17)

# ----------------------------------------------------------------- PANEL A
axA = fig.add_subplot(gs[0, 0])
h, rows = read_tsv(INTER)
I = {c: i for i, c in enumerate(h)}

def num(r, k, cast=float, default=None):
    try:
        return cast(r[I[k]])
    except (ValueError, KeyError, IndexError):
        return default

recs = []
whole = None
for r in rows:
    g = r[I["genus"]] if I["genus"] < len(r) else ""
    d = {"genus": g,
         "na": num(r, "n_amphibian", int, 0),
         "no": num(r, "n_other", int, 0),
         "obs": num(r, "obs_changes", float, 0.0),
         "nm": num(r, "null_mean", float, 0.0),
         "lo": num(r, "null_lo", float, 0.0),
         "hi": num(r, "null_hi", float, 0.0),
         "ratio": num(r, "ratio", float, None),
         "p": num(r, "p", float, None)}
    if d["ratio"] is None:
        continue
    if g == "WHOLE_TREE":
        whole = d
    else:
        recs.append(d)

recs.sort(key=lambda d: d["obs"])
n_pass = sum(1 for d in recs if d["obs"] >= MIN_TRANSITIONS)

for y, d in enumerate(recs):
    keep = d["obs"] >= MIN_TRANSITIONS
    c = AMPH if keep else GREY
    if d["nm"]:
        axA.plot([d["lo"] / d["nm"], d["hi"] / d["nm"]], [y, y],
                 color=GREY, lw=5, alpha=0.35, solid_capstyle="butt", zorder=1)
    axA.scatter([d["ratio"]], [y], s=66, color=c, zorder=3,
                edgecolor="white", linewidth=0.8)
    axA.text(-0.035, y, d["genus"] + "  (" + str(d["na"]) + "/" + str(d["no"]) + ")",
             ha="right", va="center", fontsize=8.4,
             color=DARK if keep else "#6d6d6d",
             fontweight="bold" if keep else "normal")
    axA.text(1.40, y, str(int(d["obs"])), ha="center", va="center", fontsize=8.4,
             color=DARK if keep else "#6d6d6d",
             fontweight="bold" if keep else "normal")

axA.axvline(1.0, color=DARK, lw=1.0, ls="--", zorder=2)
axA.set_xlim(-0.02, 1.50)
axA.set_ylim(-1.7, len(recs) - 0.2)
axA.set_yticks([])
axA.set_xlabel("observed / null transitions", fontsize=9)
axA.set_title("A   Host origin gives multiple independent\n"
              "transitions in " + str(n_pass) + " of " + str(len(recs)) + " genera",
              fontsize=10, loc="left", pad=9)
axA.spines["top"].set_visible(False)
axA.spines["right"].set_visible(False)
axA.spines["left"].set_visible(False)
axA.tick_params(axis="both", labelsize=8)
axA.text(1.40, len(recs) - 0.35, "transitions", ha="center", va="bottom",
         fontsize=7.6, color=DARK)
axA.text(1.02, -0.75, "no clustering", fontsize=7.4, color="#6d6d6d")
if whole is not None:
    axA.axhline(-0.75, color="#dddddd", lw=0.9)
    axA.scatter([whole["ratio"]], [-1.25], s=46, color="#555555",
                marker="D", zorder=3, edgecolor="white", linewidth=0.7)
    axA.text(-0.035, -1.25, "whole tree  (" + str(whole["na"]) + "/"
             + str(whole["no"]) + ")", ha="right", va="center",
             fontsize=8.0, color="#555555", style="italic")
    axA.text(1.40, -1.25, str(int(whole["obs"])), ha="center", va="center",
             fontsize=8.0, color="#555555")

# ----------------------------------------------------------------- PANEL B
axB = fig.add_subplot(gs[0, 1])
h, rows = read_tsv(HAPPI)
J = {c: i for i, c in enumerate(h)}
pts = []
for r in rows:
    try:
        pts.append((r[J["og"]], float(r[J["diff"]]), float(r[J["q"]])))
    except (ValueError, IndexError):
        continue

QCUT = 0.05
FLOOR = 1e-16
def ny(q):
    return -math.log10(max(q, FLOOR))

ns = [x for x in pts if x[2] >= QCUT]
sa = [x for x in pts if x[2] < QCUT and x[1] > 0]
sr = [x for x in pts if x[2] < QCUT and x[1] < 0]

axB.scatter([d for _, d, _ in ns], [ny(q) for _, _, q in ns],
            s=7, color="#d8d8d8", linewidth=0, zorder=1)
axB.scatter([d for _, d, _ in sr], [ny(q) for _, _, q in sr],
            s=11, color=REF, linewidth=0, alpha=0.75, zorder=2)
axB.scatter([d for _, d, _ in sa], [ny(q) for _, _, q in sa],
            s=11, color=AMPH, linewidth=0, alpha=0.75, zorder=2)
axB.axhline(ny(QCUT), color=DARK, lw=0.9, ls="--", zorder=3)
axB.axvline(0, color="#bbbbbb", lw=0.8, zorder=0)
axB.set_xlabel("prevalence difference (amphibian minus reference)", fontsize=9)
axB.set_ylabel("-log10 BH-adjusted p", fontsize=9)
axB.set_title("B   " + str(len(sa) + len(sr)) + " of " + str(len(pts))
              + " orthologous groups differ\n" + str(len(sa))
              + " amphibian-higher, " + str(len(sr)) + " reference-higher",
              fontsize=10, loc="left", pad=9)
axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)
axB.tick_params(axis="both", labelsize=8)
axB.text(0.97, 0.06, "amphibian-higher", transform=axB.transAxes,
         ha="right", fontsize=8.2, color=AMPH, fontweight="bold")
axB.text(0.03, 0.06, "reference-higher", transform=axB.transAxes,
         ha="left", fontsize=8.2, color=REF, fontweight="bold")

# ----------------------------------------------------------------- PANEL C
axC = fig.add_subplot(gs[0, 2])
with open(HMM, newline="") as fh:
    hrows = [r for r in csv.reader(fh, delimiter="\t") if r and len(r) >= 2]

def arm(n):
    if "__amphibian__" in n:
        return "amphibian"
    if "__endotherm__" in n or "__endo__" in n:
        return "endotherm"
    return "reference"

def gen(n):
    return n.split("__")[0] if "__" in n else "unknown"

GROUPS = [("amphibian", "Anaerotruncus"), ("reference", "Anaerotruncus"),
          ("amphibian", "UBA866"), ("reference", "UBA866")]
GLAB = ["Anaerotruncus\namphibian", "Anaerotruncus\nreference",
        "UBA866\namphibian", "UBA866\nreference"]
if not DROP_ENDO:
    GROUPS.append(("endotherm", "Anaerotruncus"))
    GLAB.append("Anaerotruncus\nendotherm-\nlabelled")

genomes = set(r[0] for r in hrows)
if os.path.exists(Q):
    qh, qr = read_tsv(Q)
    qi = qh.index("Name") if "Name" in qh else 0
    genomes |= set(r[qi].strip() for r in qr if qi < len(r))

denom = {}
for g in genomes:
    denom[(arm(g), gen(g))] = denom.get((arm(g), gen(g)), 0) + 1

carr = {}
for r in hrows:
    carr.setdefault(r[1].strip(), set()).add(r[0].strip())

ROWS = [("GATase", "GATase domain"), ("Gate", "Gate domain"),
        ("CTP_synth_N", "CTP_synth_N (pyrG)"), ("NDK", "NDK"),
        ("Nuc_H_symport", "Nuc_H_symport (NupG)")]
ROWS = [(m, l) for m, l in ROWS if m in carr]

for ri, (model, label) in enumerate(ROWS):
    y = len(ROWS) - 1 - ri
    for ci, key in enumerate(GROUPS):
        n = denom.get(key, 0)
        k = sum(1 for g in carr[model] if (arm(g), gen(g)) == key)
        f = (k / n) if n else 0.0
        s = 1.0 - 0.82 * f
        axC.add_patch(Rectangle((ci, y), 1, 1, facecolor=(s, s, s),
                      edgecolor="white", linewidth=1.6))
        axC.text(ci + 0.5, y + 0.5, str(k) + "/" + str(n), ha="center",
                 va="center", fontsize=8.2,
                 color="white" if f > 0.55 else DARK,
                 fontweight="bold" if f > 0.55 else "normal")
    axC.text(-0.12, y + 0.5, label, ha="right", va="center", fontsize=8.4)

axC.set_xlim(0, len(GROUPS))
axC.set_ylim(0, len(ROWS))
axC.set_xticks([i + 0.5 for i in range(len(GROUPS))])
axC.set_xticklabels(GLAB, fontsize=6.9, linespacing=1.35)
axC.set_yticks([])
for s in axC.spines.values():
    s.set_visible(False)
axC.tick_params(axis="x", length=0, pad=5)
axC.set_title("C   pyrG absence is arm-structured;\n"
              "NDK misses are not", fontsize=10, loc="left", pad=9)

fig.savefig(STEM + ".pdf", bbox_inches="tight")
fig.savefig(STEM + ".png", dpi=300, bbox_inches="tight")
print("wrote " + STEM + ".pdf and .png")

print()
print("=" * 78)
print("CAPTION REQUIREMENTS IMPLIED BY THIS FIGURE")
print("=" * 78)
print("PANEL A. Slatkin-Maddison test on the 1,652-tip joint tree, Fitch")
print("  parsimony, 499 label permutations. Grey bar = null 95% interval")
print("  scaled by the null mean. Brackets are amphibian/other genome counts.")
print("  Right column is the observed number of independent transitions.")
print("  GATE: >= " + str(MIN_TRANSITIONS) + " transitions. " + str(n_pass)
      + " of " + str(len(recs)) + " genera pass.")
print("  MUST STATE: the gate is the transition COUNT, not the ratio. Genera")
print("  with one transition are single monophyletic blocks in which host")
print("  origin and lineage are NOT separable.")
print("  MUST STATE: three small genera (Harryflintia 8/4, JAAYCI01 5/5,")
print("  Massiliimalia 3/7) are NOT significantly clustered, but their null")
print("  intervals span 2-5 transitions, so the test has no power there.")
print("  Non-significance is absence of evidence, not interleaving.")
print("  MUST STATE: the whole-tree row is amphibian vs all other tips and is")
print("  NOT a genus. It is set apart below the axis break.")
print()
print("PANEL B. happi (Trinh, Clausen and Willis 2023), asymptotic LRT, BH.")
print("  n = " + str(len(pts)) + " groups, " + str(len(sa) + len(sr))
      + " at q < 0.05.")
print("  CHECK: CHATINDEX records 3,345 groups. This file has " + str(len(pts))
      + ".")
print("  Resolve before quoting a denominator in text.")
print("  MUST STATE: two permutation controls (seeds 101, 202) returned 0 of")
print("  3,278. happi loses power when genome quality correlates with group,")
print("  which it does here (amphibian 87.6% vs reference 94.0% completeness),")
print("  so null results are weak evidence of no difference.")
print("  MUST STATE: a within-genus direction filter retains 422 groups, and")
print("  422 is large; aggregate coherent, individual entries UNCHECKED.")
print()
print("PANEL C. Direct Pfam HMM scan, job 27249500,")
print("  work/focal_genus_pangenome/pyrg_control/hmm/hmm_genome_model.tsv")
print("  MUST STATE: GATase and Gate at 125/125 are the positive control.")
print("  MUST STATE: NDK is present but too divergent for a profile threshold")
print("  (DIAMOND recovers it at 31-35% identity, 92-94% query coverage), so")
print("  its low counts are a sensitivity limit, not absence.")
print("  MUST STATE: CTP_synth_N is 24 of 125 by HMM and 25 by DIAMOND.")
print("  MUST STATE: NupG is GENUS-structured, not arm-structured (10/46 and")
print("  0/16 in Anaerotruncus vs 47/52 and 10/10 in UBA866). No")
print("  import-dependence claim is made.")
if not DROP_ENDO:
    print("  MUST STATE: the endotherm-labelled column is ONE genome left from")
    print("  a retired framing. Justify it or set DROP_ENDO = True.")
print()
print("WHOLE FIGURE. Anaerotruncus and UBA866/Paludihabitans only, 98")
print("  amphibian and 26 reference units" + (" plus 1 endotherm-labelled"
      if not DROP_ENDO else "") + ", one genome per species cluster.")
print("  CheckM2 1.1.0 run uniformly on all genomes in this set.")
print()
print("FIG_GENE_CONTENT_V2_20260806_COMPLETE")
# FIG_GENE_CONTENT_V2_20260806_COMPLETE
