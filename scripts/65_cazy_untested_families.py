# Families too rare to test in the focal set are scanned across the genus-unassigned clades.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/cazy_untested_families.py
# Output: results/cazy_untested_family_scan.tsv
import os, sys, csv
from collections import Counter, defaultdict

BASE = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
MTX  = os.path.join(BASE, "results/cazy_focal_family_matrix.tsv")
HAP  = os.path.join(BASE, "results/happi_cazy_focal.tsv")
OUT  = os.path.join(BASE, "results/cazy_untested_family_scan.tsv")

MIN_G = 10

def die(m):
    print("")
    print("!" * 72)
    print("FAILED: " + m)
    print("!" * 72)
    sys.exit(1)

for p in (MTX, HAP):
    if not os.path.exists(p):
        die("missing " + p)

rows = list(csv.DictReader(open(MTX), delimiter="\t"))
meta_cols = {"genome_id", "set", "genus", "arm", "clade", "total"}
fams = [c for c in rows[0].keys() if c not in meta_cols]
focal = [r for r in rows if r["set"] == "focal125" and r["arm"] in ("amphibian", "reference")]
unass = [r for r in rows if r["set"] == "unassigned_clade"]
print("=" * 72)
print("STEP 1  SETS")
print("=" * 72)
print("  families: %d" % len(fams))
print("  focal genomes in the test: %d  (endotherm-labelled excluded)" % len(focal))
print("  unassigned genomes: %d" % len(unass))
clades = defaultdict(list)
for r in unass:
    clades[r["clade"]].append(r)
print("  clade sizes: %s" % {k: len(v) for k, v in sorted(clades.items())})

tested = set()
sig = {}
for r in csv.DictReader(open(HAP), delimiter="\t"):
    tested.add(r["family"])
    if r["q"] not in ("", "NA") and float(r["q"]) < 0.05:
        sig[r["family"]] = float(r["diff"])
print("  tested: %d   significant: %d" % (len(tested), len(sig)))

print("")
print("=" * 72)
print("STEP 2  WHY EACH FAMILY WAS DROPPED")
print("=" * 72)
nf = len(focal)
pf = {f: sum(1 for r in focal if int(r[f]) > 0) for f in fams}
pu = {f: sum(1 for r in unass if int(r[f]) > 0) for f in fams}
untested = [f for f in fams if f not in tested]
too_rare = [f for f in untested if pf[f] < MIN_G]
too_common = [f for f in untested if pf[f] > nf - MIN_G]
print("  untested: %d" % len(untested))
print("    too rare in the focal set (<%d of %d): %d" % (MIN_G, nf, len(too_rare)))
print("    near-universal (>%d of %d): %d" % (nf - MIN_G, nf, len(too_common)))
print("  near-universal families: %s" % sorted(too_common))

print("")
print("=" * 72)
print("STEP 3  RARE IN FOCAL, COMMON IN THE UNASSIGNED CLADES")
print("=" * 72)
nu = len(unass)
cand = []
for f in too_rare:
    prev_f = pf[f] / float(nf)
    prev_u = pu[f] / float(nu)
    if pu[f] >= 5 and prev_u > prev_f:
        cand.append((f, pf[f], prev_f, pu[f], prev_u, prev_u - prev_f))
cand.sort(key=lambda t: -t[5])
print("  families in >=5 of %d unassigned genomes and more prevalent there" % nu)
print("  than in the focal set: %d" % len(cand))
if not cand:
    print("  NONE. The untested families are rare in both sets, so nothing is")
    print("  hidden in the white columns.")
else:
    print("")
    print("  %-24s %10s %10s %10s %10s" % ("family", "focal n", "focal %", "unass n", "unass %"))
    for f, a, pa, b, pb, d in cand[:30]:
        print("  %-24s %10d %9.1f%% %10d %9.1f%%" % (f, a, 100 * pa, b, 100 * pb))
    if len(cand) > 30:
        print("  ... and %d more" % (len(cand) - 30))

print("")
print("=" * 72)
print("STEP 4  ARE THEY CLADE-SPECIFIC OR SPREAD?")
print("=" * 72)
if cand:
    print("  A family in ONE clade is one observation. Spread across clades is")
    print("  the only thing that would be interesting, and even that is post hoc.")
    print("")
    for f, a, pa, b, pb, d in cand[:15]:
        hits = Counter()
        for c, rs in clades.items():
            n = sum(1 for r in rs if int(r[f]) > 0)
            if n:
                hits[c] = "%d/%d" % (n, len(rs))
        print("  %-24s in %d of %d clades: %s"
              % (f, len(hits), len(clades), dict(hits)))

print("")
print("=" * 72)
print("STEP 5  THE REVERSE, common in focal and absent from the clades")
print("=" * 72)
rev = []
for f in fams:
    prev_f = pf[f] / float(nf)
    prev_u = pu[f] / float(nu)
    if prev_f >= 0.30 and prev_u <= 0.05:
        rev.append((f, pf[f], prev_f, pu[f], prev_u))
rev.sort(key=lambda t: -t[2])
print("  families in >=30%% of focal genomes but <=5%% of the unassigned: %d" % len(rev))
for f, a, pa, b, pb in rev[:20]:
    mark = " [SIGNIFICANT]" if f in sig else ("" if f in tested else " [untested]")
    print("  %-24s focal %9.1f%%   unassigned %9.1f%%%s" % (f, 100 * pa, 100 * pb, mark))

if os.path.exists(OUT):
    print("")
    print("  NOT overwriting existing " + OUT)
else:
    with open(OUT, "w") as fh:
        fh.write("family\tclass\tfocal_n\tfocal_prev\tunassigned_n\tunassigned_prev\t"
                 "tested\tsignificant\tdrop_reason\n")
        for f in fams:
            reason = ""
            if f not in tested:
                reason = "too_rare" if pf[f] < MIN_G else "near_universal"
            cl = next((p for p in ("GH", "GT", "PL", "CE", "AA", "CBM") if f.startswith(p)), "other")
            fh.write("%s\t%s\t%d\t%.4f\t%d\t%.4f\t%d\t%d\t%s\n"
                     % (f, cl, pf[f], pf[f] / float(nf), pu[f], pu[f] / float(nu),
                        int(f in tested), int(f in sig), reason))
    print("")
    print("  wrote " + OUT)

print("")
print("  THIS IS A POST HOC SCAN, NOT A TEST. The unassigned set has no matched")
print("  comparison group, n is 29 against 125, and clade sizes run 2 to 11.")
print("  Nothing here gets a p value and nothing here is a result on its own.")
print("")
print("CAZY_UNTESTED_FAMILIES_V1_20260806 COMPLETE")
# CAZY_UNTESTED_FAMILIES_V1_20260806
