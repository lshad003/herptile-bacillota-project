# Contiguity is measured from the staged fastas so that all arms are measured the same way.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/assembly_quality_arms.py
# Output: results/assembly_quality_arms.tsv
# ASSEMBLY_QUALITY_ARMS_V1_20260806
# Computes N50, contig count and assembly length for every genome in both
# antiSMASH runs, from the staged fastas, so all three arms are measured the
# same way. No N50 exists in any manifest and contig counts exist only for
# the 718 wild reps (gunc_audit_by_sgb.tsv), so mixing sources would itself
# be a confound.
#
# WHY: the 2026-08-06 literature review ranked as NOT SUPPORTED any
# arithmetic decomposition of a between-arm complete-BGC gap that treats
# total regions as a fragmentation-free baseline. Total region count is
# inflated by BGC splitting across contigs and deflated by antiSMASH's
# detection floor, so it is not a stable denominator. The earlier claim that
# fragmentation explains about 80% of the gap is withdrawn.
#
# The defensible move is to restrict inference to arms that are matched on
# assembly quality and to treat GTDB as taxonomic background only. This
# script TESTS whether the two MAG arms are actually matched, rather than
# assuming it from their similar contig-edge rates (85.8% vs 84.4%).
#
# Standardized mean difference: |mean_a - mean_b| / pooled SD.
# Convention in matching work is |SMD| < 0.1 well balanced, < 0.25 tolerable.
# Also reports a two-sample Kolmogorov-Smirnov test on the distributions,
# since SMD only compares means.
#
# NOTE A SECOND CONFOUND THE REVIEW RAISED: the UHM arm is dRep
# representatives, and dRep selects on completeness, contamination and N50,
# so that arm is quality-enriched BY CONSTRUCTION. The EHI newt arm is not
# dereplicated. This pushes opposite to fragmentation. Report both.

import os, sys, gzip
import numpy as np

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
DIRS = [(B + "/work/bgc/genomes", "amphibian_run"),
        (B + "/work/bgc_refs/genomes", "reference_run")]
GEN_A = B + "/results/bgc_per_genome.tsv"
GEN_R = B + "/results/bgc_per_genome_refs.tsv"
OUT = B + "/results/assembly_quality_arms.tsv"

if os.path.exists(OUT):
    raise SystemExit("REFUSING TO OVERWRITE %s, move it first" % OUT)

try:
    from scipy.stats import ks_2samp
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False
print("scipy available:", HAVE_SCIPY)


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


def contig_lengths(path):
    lens = []
    n = 0
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                lens.append(0)
            elif lens:
                lens[-1] += len(line.strip())
    return [x for x in lens if x > 0]


def n50(lens):
    if not lens:
        return 0
    s = sorted(lens, reverse=True)
    half = sum(s) / 2.0
    run = 0
    for L in s:
        run += L
        if run >= half:
            return L
    return s[-1]


qual = {}
for d, tag in DIRS:
    if not os.path.isdir(d):
        sys.exit("MISSING: %s" % d)
    files = [e.name for e in os.scandir(d) if e.is_file()]
    print("%s: %d files" % (tag, len(files)))
    for i, fn in enumerate(sorted(files)):
        gid = fn
        for suf in (".fna.gz", ".fa.gz", ".fna", ".fa"):
            if gid.endswith(suf):
                gid = gid[: -len(suf)]
                break
        lens = contig_lengths(os.path.join(d, fn))
        qual[gid] = dict(n_contigs=len(lens), total_bp=sum(lens),
                         n50=n50(lens),
                         longest=max(lens) if lens else 0)
        if (i + 1) % 500 == 0:
            print("   %d done" % (i + 1))
print("genomes measured: %d" % len(qual))

gen_a = read_tsv(GEN_A)
gen_r = read_tsv(GEN_R)

LABEL = {("herptile", "amphibian"): "UHM amphibian",
         ("herptile", "reptile"): "UHM reptile",
         ("ehi_amphibian", "amphibian"): "EHI newt"}

rows = []
missing = 0
for r in gen_a:
    if r.get("family", "") != "Ruminococcaceae":
        continue
    arm = LABEL.get((r.get("arm", ""), r.get("host", "")))
    if arm is None:
        continue
    q = qual.get(r["genome"])
    if q is None:
        missing += 1
        continue
    rows.append(dict(genome=r["genome"], arm=arm, genus=r.get("genus", ""),
                     n_regions=int(r["n_regions"]), n_complete=int(r["n_complete"]),
                     n_edge=int(r["n_edge"]), **q))
for r in gen_r:
    q = qual.get(r["genome"])
    if q is None:
        missing += 1
        continue
    rows.append(dict(genome=r["genome"], arm="GTDB reference",
                     genus=r.get("genus", ""),
                     n_regions=int(r["n_regions"]), n_complete=int(r["n_complete"]),
                     n_edge=int(r["n_edge"]), **q))
print("rows joined: %d, genomes with no fasta measured: %d" % (len(rows), missing))

ARMS = ["UHM amphibian", "EHI newt", "GTDB reference"]
byarm = {a: [r for r in rows if r["arm"] == a] for a in ARMS}

print()
print("=" * 78)
print("ASSEMBLY QUALITY BY ARM, Ruminococcaceae only")
print("=" * 78)
print("  %-16s %6s %12s %12s %12s %12s"
      % ("arm", "n", "med N50", "med contigs", "med Mb", "med longest"))
for a in ARMS:
    v = byarm[a]
    if not v:
        continue
    print("  %-16s %6d %12.0f %12.0f %12.2f %12.0f"
          % (a, len(v),
             np.median([r["n50"] for r in v]),
             np.median([r["n_contigs"] for r in v]),
             np.median([r["total_bp"] for r in v]) / 1e6,
             np.median([r["longest"] for r in v])))

print()
print("  mean values, for the SMD calculation:")
print("  %-16s %12s %12s %12s" % ("arm", "mean N50", "mean contigs", "mean Mb"))
for a in ARMS:
    v = byarm[a]
    if not v:
        continue
    print("  %-16s %12.0f %12.1f %12.2f"
          % (a, np.mean([r["n50"] for r in v]),
             np.mean([r["n_contigs"] for r in v]),
             np.mean([r["total_bp"] for r in v]) / 1e6))


def smd(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    sp = np.sqrt((x.var(ddof=1) + y.var(ddof=1)) / 2.0)
    return abs(x.mean() - y.mean()) / sp if sp > 0 else 0.0


print()
print("=" * 78)
print("MATCHING ADEQUACY: are the two MAG arms actually comparable?")
print("=" * 78)
PAIRS = [("UHM amphibian", "EHI newt"),
         ("UHM amphibian", "GTDB reference"),
         ("EHI newt", "GTDB reference")]
METRICS = [("n50", "N50"), ("n_contigs", "contig count"),
           ("total_bp", "assembly length")]
for a, b in PAIRS:
    va, vb = byarm[a], byarm[b]
    if not va or not vb:
        continue
    print()
    print("  %s  vs  %s   (n=%d vs %d)" % (a, b, len(va), len(vb)))
    for key, lab in METRICS:
        x = [r[key] for r in va]
        y = [r[key] for r in vb]
        s = smd(x, y)
        verdict = ("well balanced" if s < 0.1 else
                   "tolerable" if s < 0.25 else "NOT BALANCED")
        line = "    %-18s SMD %6.3f   %s" % (lab, s, verdict)
        if HAVE_SCIPY:
            st, p = ks_2samp(x, y)
            line += "   KS D=%.3f p=%.2ganything" % (st, p)
            line = line.replace("anything", "")
        print(line)

print()
print("  |SMD| < 0.10 well balanced, < 0.25 tolerable, above that not matched.")
print("  SMD compares MEANS only; the KS test compares whole distributions.")
print("  If the two MAG arms are balanced, a between-arm comparison of")
print("  complete BGCs is defensible for THOSE TWO ONLY. GTDB stays as")
print("  taxonomic background regardless of what its SMD shows, because it is")
print("  MAG-heavy with study-of-origin batch effects (see R7).")

print()
print("=" * 78)
print("DEREPLICATION CHECK")
print("=" * 78)
print("  The UHM arm is dRep representatives; dRep scores on completeness,")
print("  contamination and 0.5*log10(N50), so it is quality-enriched BY")
print("  CONSTRUCTION. The EHI newt arm is NOT dereplicated. If UHM N50 is")
print("  HIGHER than EHI newt, that selection is visible and pushes opposite")
print("  to fragmentation. Report it either way.")
for a in ("UHM amphibian", "EHI newt"):
    v = byarm[a]
    if v:
        print("  %-16s median N50 %8.0f  median contigs %6.0f"
              % (a, np.median([r["n50"] for r in v]),
                 np.median([r["n_contigs"] for r in v])))

print()
print("=" * 78)
print("COMPLETE BGCs vs CONTIGUITY, WITHIN EACH ARM")
print("=" * 78)
print("  N50 quartile within arm, mean complete BGCs per genome:")
print("  %-16s %8s %8s %8s %8s" % ("arm", "Q1 low", "Q2", "Q3", "Q4 high"))
for a in ARMS:
    v = byarm[a]
    if len(v) < 20:
        continue
    ns = np.array([r["n50"] for r in v], float)
    cs = np.array([r["n_complete"] for r in v], float)
    qs = np.quantile(ns, [0.25, 0.5, 0.75])
    bins = [cs[ns <= qs[0]], cs[(ns > qs[0]) & (ns <= qs[1])],
            cs[(ns > qs[1]) & (ns <= qs[2])], cs[ns > qs[2]]]
    print("  %-16s %8.2f %8.2f %8.2f %8.2f"
          % (a, bins[0].mean(), bins[1].mean(), bins[2].mean(), bins[3].mean()))
print()
print("  A rising trend WITHIN an arm is direct evidence that contiguity")
print("  drives the complete-BGC count on this data, measured rather than")
print("  assumed. It does NOT license decomposing the between-arm gap: the")
print("  literature has no calibration converting a contiguity difference")
print("  into an expected complete-BGC difference.")

with open(OUT, "w") as f:
    f.write("genome\tarm\tgenus\tn_contigs\ttotal_bp\tn50\tlongest_contig\t"
            "n_regions\tn_complete\tn_edge\n")
    for r in sorted(rows, key=lambda z: (z["arm"], z["genome"])):
        f.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
                % (r["genome"], r["arm"], r["genus"], r["n_contigs"],
                   r["total_bp"], r["n50"], r["longest"],
                   r["n_regions"], r["n_complete"], r["n_edge"]))
print()
print("wrote %s" % OUT)
print("ASSEMBLY_QUALITY_ARMS_V1_20260806_COMPLETE")
