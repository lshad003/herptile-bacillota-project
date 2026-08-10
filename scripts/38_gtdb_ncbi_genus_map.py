# GTDB genera are mapped onto NCBI genus names, requiring a minimum share rather than any nonzero count.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/gtdb_ncbi_genus_map_v2.py
# Output: results/gtdb_ncbi_genus_map_v2.tsv
import os, csv
from collections import defaultdict, Counter

# GTDB_NCBI_GENUS_MAP V2. V1 counted any nonzero NCBI genus as a valid
# correspondence, so one mislabelled submission made a genus ambiguous.
# That put Angelakisella (293 vs three singletons) in the ambiguous bucket
# and gave NCBI Ruminococcus 141 GTDB genera including Akkermansia.
# V2 requires a minimum SHARE. Denominator is genomes carrying a named NCBI
# genus, because genomes with no NCBI genus say nothing about correspondence.

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
META = B + "/data/gtdb/bac120_metadata.tsv"
SGB = B + "/data/sgb_manifest.tsv"
OUT = B + "/results/gtdb_ncbi_genus_map_v2.tsv"

THRESHOLDS = (0.01, 0.05, 0.10)
PRIMARY = 0.05
MIN_N_CONFIDENT = 10

if os.path.exists(OUT):
    raise SystemExit("REFUSING TO OVERWRITE: " + OUT)
if not os.path.exists(META):
    raise SystemExit("MISSING: " + META)

def rank(tax, pre):
    for part in tax.split(";"):
        part = part.strip()
        if part.startswith(pre):
            return part[len(pre):].strip()
    return ""

with open(META, errors="replace") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    try:
        gi = hdr.index("gtdb_taxonomy")
        ni = hdr.index("ncbi_taxonomy")
    except ValueError:
        raise SystemExit("columns not found in header of %d fields" % len(hdr))

    pair = defaultdict(Counter)
    total = Counter()
    gfam = {}
    n = 0
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(gi, ni):
            continue
        gg = rank(p[gi], "g__")
        if not gg:
            continue
        n += 1
        total[gg] += 1
        gfam[gg] = rank(p[gi], "f__")
        ng = rank(p[ni], "g__")
        if ng:
            pair[gg][ng] += 1

print("streamed %s" % META)
print("  genomes with a GTDB genus: %d" % n)
print("  distinct GTDB genera:      %d" % len(total))
print()

# ------------------------------------------------------------- CATALOG
with open(SGB, newline="") as fh:
    r = [x for x in csv.reader(fh, delimiter="\t") if x]
J = {c: i for i, c in enumerate(r[0])}
wild = set()
for p in r[1:]:
    g = p[J["genus"]].strip()
    if g and g.upper() != "UNASSIGNED" and p[J["has_wild"]].strip() == "yes":
        wild.add(g)
print("  wild catalog genera: %d" % len(wild))
print()

def classify(gg, thr):
    c = pair.get(gg)
    if not c:
        return "no_ncbi_genus", []
    named = sum(c.values())
    keep = [(k, v) for k, v in c.most_common() if (v / named) >= thr]
    if not keep:
        return "no_ncbi_genus", []
    if len(keep) == 1:
        return "one_to_one", keep
    return "one_to_many", keep

# ------------------------------------------------------ THRESHOLD SWEEP
print("=" * 78)
print("THRESHOLD SENSITIVITY, %d WILD CATALOG GENERA" % len(wild))
print("=" * 78)
print("  %-10s %12s %12s %14s %12s" % ("min share", "one_to_one",
      "one_to_many", "no_ncbi_genus", "not_in_gtdb"))
sweep = {}
for thr in THRESHOLDS:
    cnt = Counter()
    for gg in wild:
        if gg not in total:
            cnt["absent"] += 1
            continue
        m, _ = classify(gg, thr)
        cnt[m] += 1
    sweep[thr] = cnt
    print("  %-10s %12d %12d %14d %12d"
          % ("%.0f%%" % (100 * thr), cnt["one_to_one"], cnt["one_to_many"],
             cnt["no_ncbi_genus"], cnt["absent"]))
print()
print("  TESTABLE BY NAME = one_to_one. At %.0f%% that is %d of %d."
      % (100 * PRIMARY, sweep[PRIMARY]["one_to_one"], len(wild)))
print("  V1 reported 23 with a nonzero-count rule.")
print()

# ------------------------------------------------- LOW-COUNT WARNING
low = [gg for gg in wild if 0 < total.get(gg, 0) < MIN_N_CONFIDENT]
print("  WILD GENERA WITH FEWER THAN %d REFERENCE GENOMES: %d"
      % (MIN_N_CONFIDENT, len(low)))
print("  A share threshold is meaningless at n=2. These are flagged")
print("  low_confidence in the output file.")
print()

# ------------------------------------------------- WHAT CHANGED
print("=" * 78)
print("GENERA THAT MOVE FROM ambiguous TO unambiguous BETWEEN 1%% AND %.0f%%"
      % (100 * PRIMARY))
print("=" * 78)
shown = 0
for gg in sorted(wild):
    if gg not in total:
        continue
    a, _ = classify(gg, 0.01)
    b, keep = classify(gg, PRIMARY)
    if a == "one_to_many" and b == "one_to_one" and shown < 20:
        c = pair[gg]
        print("  %-24s n=%-5d -> %s   (dropped: %s)"
              % (gg, total[gg], keep[0][0],
                 ", ".join("%s:%d" % (k, v) for k, v in c.most_common()[1:5])))
        shown += 1
print()

print("  STILL AMBIGUOUS AT %.0f%%, genuinely polyphyletic:" % (100 * PRIMARY))
shown = 0
for gg in sorted(wild):
    if gg not in total:
        continue
    m, keep = classify(gg, PRIMARY)
    if m == "one_to_many" and shown < 20:
        named = sum(pair[gg].values())
        print("  %-24s n=%-5d -> %s" % (gg, total[gg],
              ", ".join("%s %.0f%%" % (k, 100.0 * v / named) for k, v in keep)))
        shown += 1
print()

# ------------------------------------------------- REVERSE DIRECTION
print("=" * 78)
print("REVERSE DIRECTION AT %.0f%%. V1 gave Ruminococcus 141 GTDB genera."
      % (100 * PRIMARY))
print("=" * 78)
rev = defaultdict(list)
for gg in total:
    m, keep = classify(gg, PRIMARY)
    for k, v in keep:
        rev[k].append(gg)
for k in ("Ruminococcus", "Clostridium", "Eubacterium", "Blautia",
          "Faecalibacterium", "Oscillibacter", "Anaerotruncus",
          "Ruthenibacterium", "Flavonifractor", "Subdoligranulum",
          "Ethanoligenens", "Mageeibacillus", "Fastidiosipila"):
    gs = sorted(rev.get(k, []))
    print("  NCBI %-18s -> %3d GTDB genera%s"
          % (k, len(gs), (": " + ", ".join(gs[:6])
             + (" ..." if len(gs) > 6 else "")) if gs else ""))
print()

# ------------------------------------------------- WRITE
f = open(OUT, "w")
f.write("gtdb_genus\tgtdb_family\tin_wild_catalog\tn_reference_genomes\t"
        "n_with_ncbi_genus\tfrac_no_ncbi_name\tlow_confidence\t"
        "class_t01\tclass_t05\tclass_t10\tncbi_genera_at_t05\t"
        "dominant_ncbi\tdominant_share\n")
for gg in sorted(set(total) | wild):
    tot = total.get(gg, 0)
    named = sum(pair[gg].values()) if gg in pair else 0
    if tot == 0:
        f.write("%s\t\tyes\t0\t0\t0.000\tyes\tnot_in_gtdb\tnot_in_gtdb\t"
                "not_in_gtdb\t\t\t0.000\n" % gg)
        continue
    c01, _ = classify(gg, 0.01)
    c05, k05 = classify(gg, 0.05)
    c10, _ = classify(gg, 0.10)
    dom, doms = (k05[0][0], k05[0][1] / named) if k05 else ("", 0.0)
    f.write("%s\t%s\t%s\t%d\t%d\t%.3f\t%s\t%s\t%s\t%s\t%s\t%s\t%.3f\n"
            % (gg, gfam.get(gg, ""), "yes" if gg in wild else "no", tot, named,
               1.0 - (named / tot) if tot else 0.0,
               "yes" if tot < MIN_N_CONFIDENT else "no",
               c01, c05, c10,
               ";".join("%s:%.3f" % (k, v / named) for k, v in k05),
               dom, doms))
f.close()
print("  wrote %s" % OUT)
print()
print("  NEXT: the Bracken join uses class_t05 == one_to_one as the testable")
print("  set and dominant_ncbi as the name to match. Ambiguous and")
print("  no_ncbi_genus genera are reported as UNTESTABLE, not as failures.")
print()
print("GTDB_NCBI_GENUS_MAP_V2_20260806_COMPLETE")
# GTDB_NCBI_GENUS_MAP_V2_20260806_COMPLETE
