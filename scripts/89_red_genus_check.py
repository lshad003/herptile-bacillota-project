# Relative evolutionary divergence is compared between genomes with and without a GTDB genus assignment.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/red_genus_check_v2.py
# Output: results/red_genus_check_v2.tsv
import os, sys, csv
from collections import Counter

BASE = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
WILD = os.path.join(BASE, "results/gtdbtk_wild_sgb_r220/gtdbtk.bac120.summary.tsv")
AMPH = os.path.join(BASE, "results/gtdbtk_ehi_amphibian_r220/gtdbtk.bac120.summary.tsv")
COH  = os.path.join(BASE, "results/unassigned_clade_coherence.tsv")
AUD  = os.path.join(BASE, "results/seqcode_representative_audit_v2.tsv")
OUT  = os.path.join(BASE, "results/red_genus_check_v2.tsv")

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

def genus_of(tax):
    for f in str(tax).split(";"):
        f = f.strip()
        if f.startswith("g__"):
            return f[3:].strip()
    return ""

def qs(vals):
    v = sorted(vals)
    def q(f):
        return v[min(len(v) - 1, int(f * len(v)))]
    return q(0.05), q(0.25), q(0.50), q(0.75), q(0.95), v[0], v[-1]

print("=" * 72)
print("STEP 1  GTDB-Tk SUMMARIES, FORCED TAB DELIMITED")
print("=" * 72)
recs = {}
for path, label in ((WILD, "wild 718"), (AMPH, "EHI newt 437")):
    if not os.path.exists(path):
        print("  ABSENT: %s" % path)
        continue
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    if not rows:
        print("  EMPTY: %s" % path)
        continue
    hdr = [h for h in rows[0].keys() if h is not None]
    print("")
    print("  FILE: %s   (%s)" % (path, label))
    print("    rows: %d   fields: %d" % (len(rows), len(hdr)))
    for need in ("user_genome", "classification", "red_value"):
        if need not in hdr:
            print("    fields present: %s" % hdr)
            die("column '%s' missing after tab parse" % need)
    print("    user_genome / classification / red_value all present")
    nred = 0
    for r in rows:
        g = (r["user_genome"] or "").strip()
        if not g:
            continue
        red = fnum(r["red_value"])
        if red is not None:
            nred += 1
        recs[g] = {"arm": label,
                   "tax": (r["classification"] or "").strip(),
                   "genus": genus_of(r["classification"]),
                   "red": red,
                   "method": (r.get("classification_method") or "").strip(),
                   "msa": fnum(r.get("msa_percent"))}
    print("    rows with a numeric red_value: %d of %d (%.1f%%)"
          % (nred, len(rows), 100.0 * nred / len(rows)))
    print("    classification_method counts: %s"
          % dict(Counter((r.get("classification_method") or "blank").strip()
                         for r in rows).most_common(6)))
if not recs:
    die("nothing loaded")
print("")
print("  total genomes loaded: %d" % len(recs))
print("  example ids: %s" % list(recs)[:3])

print("")
print("=" * 72)
print("STEP 2  EMPIRICAL RED CALIBRATION FROM THIS RUN")
print("=" * 72)
named   = [v for v in recs.values() if v["genus"] and v["red"] is not None]
unnamed = [v for v in recs.values() if not v["genus"] and v["red"] is not None]
print("  WITH genus + RED   : %d" % len(named))
print("  WITHOUT genus + RED: %d" % len(unnamed))
print("  WITH genus, no RED : %d"
      % sum(1 for v in recs.values() if v["genus"] and v["red"] is None))
print("  WITHOUT genus, no RED: %d"
      % sum(1 for v in recs.values() if not v["genus"] and v["red"] is None))
if len(named) < 20:
    print("")
    print("  *** FEWER THAN 20 ASSIGNED GENOMES CARRY A RED VALUE ***")
    print("  Calibration below is unreliable. Do not build a claim on it.")

lo = None
if named:
    a = qs([v["red"] for v in named])
    lo = a[0]
    print("")
    print("  RED, genomes placed IN a named genus (n=%d):" % len(named))
    print("    min %.4f  5th %.4f  25th %.4f  MEDIAN %.4f  75th %.4f  95th %.4f  max %.4f"
          % (a[5], a[0], a[1], a[2], a[3], a[4], a[6]))
if unnamed:
    b = qs([v["red"] for v in unnamed])
    print("")
    print("  RED, genomes with NO genus assignment (n=%d):" % len(unnamed))
    print("    min %.4f  5th %.4f  25th %.4f  MEDIAN %.4f  75th %.4f  95th %.4f  max %.4f"
          % (b[5], b[0], b[1], b[2], b[3], b[4], b[6]))
if named and unnamed:
    below = sum(1 for v in unnamed if v["red"] < lo)
    print("")
    print("  Unassigned genomes with RED below the assigned 5th percentile")
    print("  (%.4f): %d of %d (%.1f%%)" % (lo, below, len(unnamed), 100.0 * below / len(unnamed)))
    print("")
    print("  LOWER RED = DEEPER placement = what a novel genus looks like.")
    print("  RED inside the assigned range = within-genus depth, i.e. a novel")
    print("  SPECIES in a genus GTDB has not named, NOT a novel genus.")

print("")
print("=" * 72)
print("STEP 3  THE UNASSIGNED CLADES")
print("=" * 72)
if not os.path.exists(COH):
    die("missing " + COH)
crows = list(csv.DictReader(open(COH), delimiter="\t"))
print("  clade rows: %d" % len(crows))
pass4 = set()
if os.path.exists(AUD):
    for r in csv.DictReader(open(AUD), delimiter="\t"):
        if (r.get("pass_four") or "0").strip() == "1":
            pass4.add((r["representative"] or "").strip())
    print("  representatives passing all four assembly criteria: %d" % len(pass4))

out = []
seen = nored = nomatch = 0
tally = Counter()
for r in crows:
    clade = (r["clade"] or "").strip()
    gs = [x.strip() for x in (r["genomes"] or "").split(";") if x.strip()]
    print("")
    print("  CLADE %s  (%s tips)" % (clade, r["n_tips"]))
    d = (r.get("dist_to_nearest_assigned") or "").strip()
    if d:
        print("    patristic distance to nearest assigned tip: %s" % d[:6])
    for g in gs:
        seen += 1
        v = recs.get(g)
        p4 = "PASS4" if g in pass4 else ""
        if v is None:
            nomatch += 1
            tally["not_in_summary"] += 1
            print("    %-42s NOT IN ANY SUMMARY FILE %s" % (g, p4))
            out.append((clade, g, "", "", "", "not_in_summary", int(g in pass4)))
            continue
        if v["red"] is None:
            nored += 1
            tally["no_red"] += 1
            print("    %-42s RED absent   method=%-20s %s"
                  % (g, v["method"] or "blank", p4))
            out.append((clade, g, "", v["method"], v["tax"], "no_red", int(g in pass4)))
            continue
        verdict = "UNTESTABLE"
        if lo is not None:
            verdict = "BELOW_assigned_5th" if v["red"] < lo else "within_assigned_range"
        tally[verdict] += 1
        print("    %-42s RED %.4f  %-20s %-22s %s"
              % (g, v["red"], v["method"] or "blank", verdict, p4))
        out.append((clade, g, "%.4f" % v["red"], v["method"], v["tax"], verdict,
                    int(g in pass4)))

print("")
print("  genomes checked %d | no RED %d | not in summary %d" % (seen, nored, nomatch))
print("  verdict tally: %s" % dict(tally))
if nored:
    print("  A genome with no RED was called by topology or ANI and is simply")
    print("  untestable this way. Report separately, never as pass or fail.")

if os.path.exists(OUT):
    print("")
    print("  NOT overwriting existing " + OUT)
else:
    with open(OUT, "w") as fh:
        fh.write("clade\tgenome\tred_value\tclassification_method\tclassification\tverdict\tpass_four\n")
        for row in out:
            fh.write("\t".join(str(x) for x in row) + "\n")
    print("")
    print("  wrote " + OUT)

print("")
print("=" * 72)
print("HOW TO READ THIS")
print("=" * 72)
print("  Calibrated against genomes classified in the SAME job. GTDB's own")
print("  published RED interval for genus is a different number and is NOT")
print("  used here; it comes from the r220 release files if you want it.")
print("  RED is computed on the GTDB reference tree during placement, NOT on")
print("  work/rep_tree/figure_tree.nwk, so this and the patristic clade result")
print("  are INDEPENDENT lines of evidence.")
print("")
print("RED_GENUS_CHECK_V2_20260806 COMPLETE")
# RED_GENUS_CHECK_V2_20260806
