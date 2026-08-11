# Product class composition is compared across arms, with classes excluded where the proportion tracks contiguity within arms.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/bgc_class_composition.py
# Output: results/bgc_class_composition.tsv
import os, csv, math
from collections import defaultdict, Counter

# BGC_CLASS_COMPOSITION V1. Complete-BGC COUNTS cannot be compared across
# arms: the arms are not matched on contiguity and converge at matched N50.
# Class COMPOSITION is a different question and is far less confounded,
# because fragmentation reduces how many clusters are complete without much
# changing which product classes are present.
# THE GATE: if a class proportion correlates with N50 WITHIN arms, any
# between-arm difference in that class is a fragmentation artifact and is
# not reportable. That test is run first and its verdict is printed per class.

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
AMPH = B + "/results/bgc_regions.tsv"
REFS = B + "/results/bgc_regions_refs.tsv"
ENDO = B + "/results/bgc_regions_endo.tsv"
QUAL = B + "/results/assembly_quality_all_arms.tsv"
OUT = B + "/results/bgc_class_composition.tsv"

MIN_REGIONS = 2
MIN_GENOMES = 20

if os.path.exists(OUT):
    raise SystemExit("REFUSING TO OVERWRITE: " + OUT)
for p in (AMPH, REFS, ENDO, QUAL):
    if not os.path.exists(p):
        raise SystemExit("MISSING: " + p)

def read_tsv(path):
    with open(path, newline="") as fh:
        r = [x for x in csv.reader(fh, delimiter="\t") if x]
    return r[0], r[1:]

# --------------------------------------------------- N50 PER GENOME
h, rows = read_tsv(QUAL)
Q = {c: i for i, c in enumerate(h)}
n50 = {}
qarm = {}
for p in rows:
    try:
        n50[p[Q["genome"]].strip()] = float(p[Q["n50"]])
        qarm[p[Q["genome"]].strip()] = p[Q["arm"]].strip()
    except (ValueError, IndexError):
        continue
print("N50 available for %d genomes" % len(n50))

# --------------------------------------------------- REGIONS, ONE FAMILY
per_genome = defaultdict(Counter)
arm_of = {}
dropped_family = 0

h, rows = read_tsv(AMPH)
A = {c: i for i, c in enumerate(h)}
for p in rows:
    if p[A["family"]].strip() != "Ruminococcaceae":
        dropped_family += 1
        continue
    g = p[A["genome"]].strip()
    a = p[A["arm"]].strip()
    arm_of[g] = "UHM amphibian" if a.startswith("uhm") or a == "herptile" \
        else ("EHI newt" if "amphib" in a else a)
    per_genome[g][p[A["product"]].strip()] += 1

print("amphibian file: dropped %d regions outside Ruminococcaceae"
      % dropped_family)

for path, label in ((REFS, "GTDB reference"), (ENDO, None)):
    h, rows = read_tsv(path)
    C = {c: i for i, c in enumerate(h)}
    for p in rows:
        g = p[C["genome"]].strip()
        a = p[C["arm"]].strip()
        if label:
            arm_of[g] = label
        else:
            arm_of[g] = "EHI mammal" if "mam" in a else "Youngblut"
        per_genome[g][p[C["product"]].strip()] += 1

# arms are authoritative from the quality file where available
for g in list(arm_of):
    if g in qarm:
        arm_of[g] = qarm[g]

ARMS = ["UHM amphibian", "EHI newt", "EHI mammal", "Youngblut",
        "GTDB reference"]
genomes = {a: [g for g in per_genome
               if arm_of.get(g) == a and sum(per_genome[g].values()) >= MIN_REGIONS
               and g in n50] for a in ARMS}
print()
print("GENOMES WITH >= %d REGIONS AND AN N50:" % MIN_REGIONS)
for a in ARMS:
    tot = sum(1 for g in per_genome if arm_of.get(g) == a)
    print("  %-18s %5d of %5d" % (a, len(genomes[a]), tot))
print()

# --------------------------------------------------- CLASSES
allc = Counter()
for g, c in per_genome.items():
    if arm_of.get(g) in ARMS:
        for k, v in c.items():
            allc[k] += v
CLASSES = [k for k, v in allc.most_common() if v >= 50]
print("PRODUCT CLASSES WITH >= 50 REGIONS ACROSS ALL ARMS: %d" % len(CLASSES))
print("  " + ", ".join(CLASSES))
print()

def prop(g, cls):
    tot = sum(per_genome[g].values())
    return per_genome[g][cls] / tot if tot else 0.0

def mean(v):
    return sum(v) / len(v) if v else float("nan")

def spearman(x, y):
    n = len(x)
    if n < 8:
        return float("nan"), float("nan")
    def rank(v):
        idx = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[idx[j + 1]] == v[idx[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[idx[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(x), rank(y)
    mx, my = mean(rx), mean(ry)
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx <= 0 or syy <= 0:
        return float("nan"), float("nan")
    rho = sxy / math.sqrt(sxx * syy)
    t = rho * math.sqrt(max(n - 2, 1) / max(1e-12, 1 - rho * rho))
    return rho, t

# --------------------------------------------------- THE GATE
print("=" * 78)
print("GATE. Does class proportion track N50 WITHIN each arm?")
print("A class that does is FRAGMENTATION-SENSITIVE and cannot be compared")
print("between arms that differ in contiguity.")
print("=" * 78)
print("  %-30s %s" % ("class", "  ".join("%-14s" % a[:13] for a in ARMS)))
gate = {}
for cls in CLASSES:
    cells, flags = [], []
    for a in ARMS:
        gs = genomes[a]
        if len(gs) < MIN_GENOMES:
            cells.append("  n/a         ")
            continue
        rho, t = spearman([n50[g] for g in gs], [prop(g, cls) for g in gs])
        strong = (rho == rho) and abs(rho) >= 0.20 and abs(t) >= 2.0
        flags.append(strong)
        cells.append("%+0.3f%s        " % (rho, "*" if strong else " ")[:14])
    gate[cls] = sum(1 for f in flags if f)
    print("  %-30s %s" % (cls, "".join(c[:16] for c in cells)))
print()
print("  * marks |rho| >= 0.20 with |t| >= 2. Count of flagged arms per class:")
for cls in CLASSES:
    print("    %-30s %d of %d arms" % (cls, gate[cls],
          sum(1 for a in ARMS if len(genomes[a]) >= MIN_GENOMES)))
print()

# --------------------------------------------------- COMPOSITION
print("=" * 78)
print("MEAN PER-GENOME PROPORTION BY CLASS AND ARM")
print("Genomes with >= %d regions only. Per-genome proportions, NOT pooled" %
      MIN_REGIONS)
print("region counts, so cluster-rich genomes cannot dominate.")
print("=" * 78)
print("  %-30s %s  %s" % ("class",
      "  ".join("%-13s" % a[:12] for a in ARMS), "gate"))
out_rows = []
for cls in CLASSES:
    cells = []
    for a in ARMS:
        gs = genomes[a]
        m = mean([prop(g, cls) for g in gs]) if len(gs) >= MIN_GENOMES \
            else float("nan")
        cells.append("%0.3f        " % m if m == m else "  n/a        ")
        out_rows.append((cls, a, len(gs), m, gate[cls]))
    tag = "OK" if gate[cls] == 0 else "FRAGMENTATION-SENSITIVE"
    print("  %-30s %s  %s" % (cls, "".join(c[:15] for c in cells), tag))
print()

print("  AMPHIBIAN ARMS vs EHI MAMMAL, classes that PASSED the gate:")
for cls in CLASSES:
    if gate[cls] != 0:
        continue
    vals = {}
    for a in ARMS:
        if len(genomes[a]) >= MIN_GENOMES:
            vals[a] = mean([prop(g, cls) for g in genomes[a]])
    if "UHM amphibian" in vals and "EHI mammal" in vals:
        d = vals["UHM amphibian"] - vals["EHI mammal"]
        print("    %-30s UHM %0.3f  newt %0.3f  mammal %0.3f  diff %+0.3f"
              % (cls, vals.get("UHM amphibian", float("nan")),
                 vals.get("EHI newt", float("nan")),
                 vals.get("EHI mammal", float("nan")), d))
print()

f = open(OUT, "w")
f.write("product_class\tarm\tn_genomes\tmean_per_genome_proportion\t"
        "n_arms_fragmentation_flagged\n")
for cls, a, n, m, gt in out_rows:
    f.write("%s\t%s\t%d\t%s\t%d\n"
            % (cls, a, n, ("%.6f" % m) if m == m else "NA", gt))
f.close()
print("  wrote %s" % OUT)
print()
print("  READ THE GATE COLUMN FIRST. Any class flagged in one or more arms")
print("  must NOT be compared between arms, however large the difference.")
print()
print("BGC_CLASS_COMPOSITION_V1_20260806_COMPLETE")
# BGC_CLASS_COMPOSITION_V1_20260806_COMPLETE
