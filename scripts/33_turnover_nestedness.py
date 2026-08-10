# Genus pools are partitioned into turnover and nestedness, with threshold sensitivity and accumulation.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/turnover_nestedness.py
# Output: results/turnover_nestedness.tsv, threshold_sensitivity.tsv, genus_accumulation.tsv
# TURNOVER_NESTEDNESS_V1_20260804
import os, random
from collections import Counter, defaultdict

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
WORK = os.path.join(ROOT, "work/pooled_drep_rum")
CDB = os.path.join(WORK, "drep_out/data_tables/Cdb.csv")
MAP = os.path.join(WORK, "genome_arms.tsv")

HERP_SUM = os.path.join(ROOT, "results/gtdbtk_wild_sgb_r220/gtdbtk.bac120.summary.tsv")
AMPH_SUM = os.path.join(ROOT, "results/gtdbtk_ehi_amphibian_r220/gtdbtk.bac120.summary.tsv")
EHI_SUM = os.path.join(ROOT, "results/gtdbtk_ehi_r220_classify/gtdbtk.bac120.summary.tsv")
YB_SUM = os.path.join(ROOT, "results/gtdbtk_youngblut_r220/gtdbtk.bac120.summary.tsv")
GTDB_TAX = "/srv/projects/db/gtdbtk/220/taxonomy/gtdb_taxonomy.tsv"

OUT_BETA = os.path.join(ROOT, "results/turnover_nestedness.tsv")
OUT_THRESH = os.path.join(ROOT, "results/threshold_sensitivity.tsv")
OUT_ACC = os.path.join(ROOT, "results/genus_accumulation.tsv")

N_BOOT = 499
random.seed(20260804)
ARMS = ["herptile", "ehi_amphibian", "ehi", "youngblut", "gtdb_ref"]
AMPH_ARMS = ["herptile", "ehi_amphibian"]
ENDO_ARMS = ["ehi", "youngblut"]


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


# ------------------------------------------------------------ genus labels
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
            acc = f[0].strip().replace("GB_", "").replace("RS_", "")
            genus[("gtdb_ref", acc)] = norm(t.get("g", ""))

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

# one unit per arm per cluster, jointly dereplicated
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

print("=" * 74)
print("UNITS (jointly dereplicated, 95% ANI, one per arm per cluster)")
print("=" * 74)
for a in ARMS:
    print("  %-14s %5d units, %3d genera" % (a, len(units.get(a, [])), len(sets[a])))

# ============================================================ BETAPART
print()
print("=" * 74)
print("TURNOVER vs NESTEDNESS (Baselga 2010 decomposition of Sorensen)")
print("=" * 74)
print("  b_sor = total dissimilarity = (b+c)/(2a+b+c)")
print("  b_sim = TURNOVER (replacement)  = min(b,c)/(a+min(b,c))")
print("  b_sne = NESTEDNESS component    = b_sor - b_sim")
print()
print("  If turnover dominates, the pools are genuinely DIFFERENT.")
print("  If nestedness dominates, one pool is an UNDERSAMPLED SUBSET of the")
print("  other, and the disjointness claim collapses into a sampling story.")
print()


def beta(sa, sb):
    a = len(sa & sb)
    b = len(sa - sb)
    c = len(sb - sa)
    sor = (b + c) / float(2 * a + b + c) if (2 * a + b + c) else 0.0
    sim = min(b, c) / float(a + min(b, c)) if (a + min(b, c)) else 0.0
    return a, b, c, sor, sim, sor - sim


PAIRS = [("herptile", "ehi"), ("herptile", "youngblut"),
         ("ehi_amphibian", "ehi"), ("ehi_amphibian", "youngblut"),
         ("herptile", "ehi_amphibian"), ("ehi", "youngblut"),
         ("herptile", "gtdb_ref"), ("ehi", "gtdb_ref")]

out = open(OUT_BETA, "w")
out.write("arm_a\tarm_b\tshared\tonly_a\tonly_b\tb_sor\tb_sim_turnover\t"
          "b_sne_nestedness\tturnover_pct\tdominant\n")
print("  %-30s %5s %5s %5s %7s %7s %7s %8s"
      % ("pair", "share", "onlyA", "onlyB", "b_sor", "b_sim", "b_sne", "turn%"))
for x, y in PAIRS:
    a, b, c, sor, sim, sne = beta(sets[x], sets[y])
    pct = 100.0 * sim / sor if sor else float("nan")
    dom = "turnover" if sim > sne else "nestedness"
    print("  %-30s %5d %5d %5d %7.3f %7.3f %7.3f %7.1f%%  %s"
          % (x + " vs " + y, a, b, c, sor, sim, sne, pct, dom))
    out.write("%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.2f\t%s\n"
              % (x, y, a, b, c, sor, sim, sne, pct, dom))
out.close()
print()
print("  wrote %s" % OUT_BETA)

# ============================================================ POSITIVE CONTROL
print()
print("=" * 74)
print("POSITIVE CONTROL: GENERA RECOVERED IN EVERY ARM")
print("=" * 74)
print("  An absence claim needs proof the pipeline CAN recover this family")
print("  from every arm. A genus present in all five is that proof.")
in_all = sorted(set.intersection(*[sets[a] for a in ARMS]),
                key=lambda g: -sum(counts[a].get(g, 0) for a in ARMS))
print()
print("  %-24s %9s %6s %6s %6s %6s" % ("genus", "herptile", "newt", "EHImam", "Yblut", "refs"))
for g in in_all:
    print("  %-24s %9d %6d %6d %6d %6d"
          % (g, counts["herptile"].get(g, 0), counts["ehi_amphibian"].get(g, 0),
             counts["ehi"].get(g, 0), counts["youngblut"].get(g, 0),
             counts["gtdb_ref"].get(g, 0)))
print()
print("  genera in all five arms: %d" % len(in_all))
if not in_all:
    print("  NONE. Without a positive control the absence claim is not supportable.")

print()
print("  BIDIRECTIONAL CHECK: genera in both endotherm arms and neither")
print("  amphibian arm (argues against one-directional recovery bias)")
endo_only = sorted([g for g in sets["ehi"] & sets["youngblut"]
                    if g not in sets["herptile"] and g not in sets["ehi_amphibian"]],
                   key=lambda g: -(counts["ehi"].get(g, 0) + counts["youngblut"].get(g, 0)))
for g in endo_only[:10]:
    print("    %-24s EHImam %3d, Yblut %3d, refs %3d"
          % (g, counts["ehi"].get(g, 0), counts["youngblut"].get(g, 0),
             counts["gtdb_ref"].get(g, 0)))
print("    total: %d genera" % len(endo_only))

# ============================================================ THRESHOLD SWEEP
print()
print("=" * 74)
print("THRESHOLD SENSITIVITY")
print("=" * 74)
print("  How fast does the 'in both amphibian arms, in neither endotherm arm'")
print("  set erode as the minimum-genomes-to-count-as-present rises?")
print("  A set that survives only at threshold 1 is a threshold artefact.")
print()
out = open(OUT_THRESH, "w")
out.write("min_genomes\tn_genera\therptile_units\tnewt_units\tgenera_list\n")
print("  %-12s %8s %10s %8s  %s" % ("min genomes", "genera", "herp units", "newt", "lost at this step"))
prev = None
for thr in range(1, 6):
    def present(a, g):
        return counts[a].get(g, 0) >= thr
    cand = [g for g in sets["herptile"] | sets["ehi_amphibian"]
            if present("herptile", g) and present("ehi_amphibian", g)
            and not present("ehi", g) and not present("youngblut", g)]
    mh = sum(counts["herptile"].get(g, 0) for g in cand)
    mn = sum(counts["ehi_amphibian"].get(g, 0) for g in cand)
    lost = "" if prev is None else ", ".join(sorted(set(prev) - set(cand))) or "none"
    print("  %-12d %8d %10d %8d  %s" % (thr, len(cand), mh, mn, lost))
    out.write("%d\t%d\t%d\t%d\t%s\n" % (thr, len(cand), mh, mn, ";".join(sorted(cand))))
    prev = cand
out.close()
print()
print("  wrote %s" % OUT_THRESH)

# ============================================================ POWER
print()
print("=" * 74)
print("DETECTION POWER: 1-(1-p)^N")
print("=" * 74)
print("  Probability of recovering at least one genome of a genus present at")
print("  prevalence p, given a catalog of N units. This is what an absence in")
print("  a small catalog is actually worth.")
print()
print("  %-8s %s" % ("p", "  ".join("%-14s" % ("%s N=%d" % (a[:9], len(units.get(a, []))))
                                    for a in ARMS)))
for p in (0.01, 0.02, 0.05, 0.10, 0.20):
    cells = []
    for a in ARMS:
        N = len(units.get(a, []))
        cells.append("%-14.3f" % (1 - (1 - p) ** N))
    print("  %-8.2f %s" % (p, "  ".join(cells)))
print()
print("  Read the youngblut column: at low prevalence a 48-unit catalog has")
print("  poor power, so absence there is weak evidence and should not by")
print("  itself disqualify a genus from the amphibian-associated set.")

# ============================================================ ACCUMULATION
print()
print("=" * 74)
print("GENUS ACCUMULATION (has each catalog's genus inventory saturated?)")
print("=" * 74)
out = open(OUT_ACC, "w")
out.write("arm\tn_sampled\tmean_genera\tlo\thi\n")
DEPTHS = [5, 10, 20, 30, 46, 60, 80, 120, 180, 218]
print("  %-8s %s" % ("n", "  ".join("%-16s" % a[:15] for a in ARMS)))
for d in DEPTHS:
    cells = []
    for a in ARMS:
        lst = units.get(a, [])
        if d > len(lst):
            cells.append("%-16s" % ".")
            continue
        vals = [len(set(random.sample(lst, d))) for _ in range(200)]
        vals.sort()
        m = sum(vals) / float(len(vals))
        lo, hi = vals[int(0.025 * (len(vals) - 1))], vals[int(0.975 * (len(vals) - 1))]
        cells.append("%-16s" % ("%.1f [%d,%d]" % (m, lo, hi)))
        out.write("%s\t%d\t%.2f\t%d\t%d\n" % (a, d, m, lo, hi))
    print("  %-8d %s" % (d, "  ".join(cells)))
out.close()
print()
print("  wrote %s" % OUT_ACC)
print("  A curve still climbing at its final depth has NOT saturated, so its")
print("  absences are sampling-limited. Compare youngblut and newt to herptile.")

# SENTINEL_END
