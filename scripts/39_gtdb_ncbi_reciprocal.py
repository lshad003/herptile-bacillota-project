# The mapping is filtered to genera reciprocal at 70 percent in both directions, since a forward share alone is not sufficient.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/gtdb_ncbi_reciprocal.py
# Output: results/gtdb_ncbi_reciprocal.tsv
import os, csv
from collections import defaultdict, Counter

# GTDB_NCBI_RECIPROCAL V1. V2 used a FORWARD share only, which is not
# sufficient: GTDB Merdicola maps one-to-one onto NCBI Clostridium, but NCBI
# Clostridium maps onto 213 GTDB genera, so a Clostridium read cannot be
# attributed to Merdicola. V2's counts were also monotonic in the threshold
# (29/39/49), so a looser rule always produced more "testable" genera.
# Here a pair is usable only if it is RECIPROCAL: the GTDB genus is mostly
# labelled N in NCBI, AND genomes labelled N in NCBI are mostly in that
# GTDB genus.

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
META = B + "/data/gtdb/bac120_metadata.tsv"
SGB = B + "/data/sgb_manifest.tsv"
OUT = B + "/results/gtdb_ncbi_reciprocal.tsv"

THRESHOLDS = (0.50, 0.70, 0.90)
PRIMARY = 0.70

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

def norm(n):
    if n.startswith("Candidatus "):
        return n[len("Candidatus "):].strip()
    return n

with open(META, errors="replace") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    gi = hdr.index("gtdb_taxonomy")
    ni = hdr.index("ncbi_taxonomy")
    pair = Counter()
    gtot = Counter()
    ntot = Counter()
    gfam = {}
    n = 0
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(gi, ni):
            continue
        gg = rank(p[gi], "g__")
        if not gg:
            continue
        gfam[gg] = rank(p[gi], "f__")
        ng = norm(rank(p[ni], "g__"))
        if not ng:
            continue
        n += 1
        pair[(gg, ng)] += 1
        gtot[gg] += 1
        ntot[ng] += 1

print("streamed %s" % META)
print("  genomes with BOTH a GTDB genus and an NCBI genus: %d" % n)
print("  GTDB genera: %d    NCBI genera: %d" % (len(gtot), len(ntot)))
print()

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

def best(gg):
    cands = [(k[1], v) for k, v in pair.items() if k[0] == gg]
    if not cands:
        return None
    return max(cands, key=lambda x: x[1])

by_g = defaultdict(list)
for (gg, ng), v in pair.items():
    by_g[gg].append((ng, v))

def evaluate(gg, thr):
    if gg not in gtot:
        return "not_in_gtdb", "", 0.0, 0.0
    ng, v = max(by_g[gg], key=lambda x: x[1])
    fwd = v / gtot[gg]
    rev = v / ntot[ng]
    if fwd < thr and rev < thr:
        return "both_directions_fail", ng, fwd, rev
    if fwd < thr:
        return "forward_ambiguous", ng, fwd, rev
    if rev < thr:
        return "ncbi_name_is_shared", ng, fwd, rev
    return "reciprocal", ng, fwd, rev

print("=" * 78)
print("RECIPROCAL CRITERION, %d WILD CATALOG GENERA" % len(wild))
print("=" * 78)
print("  %-10s %12s %18s %20s %12s" % ("min share", "reciprocal",
      "forward_ambiguous", "ncbi_name_is_shared", "both_fail"))
for thr in THRESHOLDS:
    c = Counter()
    for gg in wild:
        m, _, _, _ = evaluate(gg, thr)
        c[m] += 1
    print("  %-10s %12d %18d %20d %12d"
          % ("%.0f%%" % (100 * thr), c["reciprocal"], c["forward_ambiguous"],
             c["ncbi_name_is_shared"], c["both_directions_fail"]))
    if c["not_in_gtdb"]:
        print("    plus %d not in GTDB metadata" % c["not_in_gtdb"])
print()
print("  Unlike the forward-only rule, this is NOT monotonic in the")
print("  direction that flatters the analysis: raising the bar tightens")
print("  BOTH sides.")
print()

print("=" * 78)
print("RECIPROCALLY TESTABLE WILD GENERA AT %.0f%%" % (100 * PRIMARY))
print("=" * 78)
ok = []
for gg in sorted(wild):
    m, ng, f_, r_ = evaluate(gg, PRIMARY)
    if m == "reciprocal":
        ok.append((gg, ng, f_, r_, gtot[gg], ntot[ng]))
print("  %-24s %-24s %7s %7s %8s %8s"
      % ("GTDB genus", "NCBI genus", "fwd", "rev", "n_gtdb", "n_ncbi"))
for gg, ng, f_, r_, a, b in sorted(ok, key=lambda x: -x[4]):
    print("  %-24s %-24s %6.0f%% %6.0f%% %8d %8d"
          % (gg, ng, 100 * f_, 100 * r_, a, b))
print()
print("  TOTAL RECIPROCAL: %d of %d wild genera" % (len(ok), len(wild)))
print()

print("  REJECTED BECAUSE THE NCBI NAME IS SHARED WITH OTHER GTDB GENERA:")
sh = 0
for gg in sorted(wild):
    m, ng, f_, r_ = evaluate(gg, PRIMARY)
    if m == "ncbi_name_is_shared" and sh < 20:
        print("    %-24s -> %-20s fwd %3.0f%% but rev only %3.0f%% (%d genomes"
              " named %s across GTDB)" % (gg, ng, 100 * f_, 100 * r_,
              ntot[ng], ng))
        sh += 1
print()

f = open(OUT, "w")
f.write("gtdb_genus\tgtdb_family\tin_wild_catalog\tbest_ncbi_genus\t"
        "n_gtdb_genomes\tn_ncbi_genomes\tn_shared\tforward_share\t"
        "reverse_share\tverdict_t50\tverdict_t70\tverdict_t90\n")
for gg in sorted(set(gtot) | wild):
    if gg not in gtot:
        f.write("%s\t\tyes\t\t0\t0\t0\t0.000\t0.000\tnot_in_gtdb\t"
                "not_in_gtdb\tnot_in_gtdb\n" % gg)
        continue
    ng, v = max(by_g[gg], key=lambda x: x[1])
    f.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%.3f\t%.3f\t%s\t%s\t%s\n"
            % (gg, gfam.get(gg, ""), "yes" if gg in wild else "no", ng,
               gtot[gg], ntot[ng], v, v / gtot[gg], v / ntot[ng],
               evaluate(gg, 0.50)[0], evaluate(gg, 0.70)[0],
               evaluate(gg, 0.90)[0]))
f.close()
print("  wrote %s" % OUT)
print()
print("GTDB_NCBI_RECIPROCAL_V1_20260806_COMPLETE")
# GTDB_NCBI_RECIPROCAL_V1_20260806_COMPLETE
