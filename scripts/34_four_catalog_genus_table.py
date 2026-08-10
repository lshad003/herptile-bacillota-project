# Genera are assigned to recovery blocks across four catalogs, with richness-matched Jaccard.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/four_catalog_genus_table.py
# Output: results/four_catalog_genus_table.tsv, four_catalog_jaccard.tsv, amphibian_only_genus_sources.tsv
# FOUR_CATALOG_GENUS_TABLE_V2_20260806
# V2 changes the BLOCK_LABEL strings ONLY. All computation is identical to V1.
#
# V1 labelled block 2 "NOT recovered from any endotherm gut". That is the
# claim R3 forbids. block_of() tests ONLY EHI mammals and Youngblut, and
# those genera have 62 GTDB reference genomes of which 42 are gut or faecal
# sources including 17 ruminant GI tract, 4 chicken caecum, 4 mouse gut,
# 2 human and 1 rat cecum (results/amphibian_only_genus_sources.tsv).
# This is a DIFFERENTIAL RECOVERY claim, not a compositional absence claim.
# Same error class as the Figure 3A label fixed in fig3_cross_catalog.py V7.
#
# NOTE ON ARMS: "ehi" here is the 2,481 EHI MAMMAL/BIRD genomes only
# (gtdbtk_ehi_r220_classify). EHI newts are NOT an arm in this script.
# The newt comparison lives in results/amphibian_genus_replication.tsv,
# which cross-tabulates in_both_amphibian against in_any_endotherm and
# gives 8 / 6 / 25 / 98. Both block schemes are correct; they answer
# different questions. Do not treat them as competing versions.

import os, random
from collections import Counter, defaultdict

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
WORK = os.path.join(ROOT, "work/pooled_drep_rum")
CDB = os.path.join(WORK, "drep_out/data_tables/Cdb.csv")
WDB = os.path.join(WORK, "drep_out/data_tables/Wdb.csv")
MAP = os.path.join(WORK, "genome_arms.tsv")

SGB_MANIFEST = os.path.join(ROOT, "data/sgb_manifest.tsv")
HERP_SUM = os.path.join(ROOT, "results/gtdbtk_wild_sgb_r220/gtdbtk.bac120.summary.tsv")
EHI_SUM = os.path.join(ROOT, "results/gtdbtk_ehi_r220_classify/gtdbtk.bac120.summary.tsv")
YB_SUM = os.path.join(ROOT, "results/gtdbtk_youngblut_r220/gtdbtk.bac120.summary.tsv")
GTDB_TAX = "/srv/projects/db/gtdbtk/220/taxonomy/gtdb_taxonomy.tsv"

OUT_TABLE = os.path.join(ROOT, "results/four_catalog_genus_table.tsv")
OUT_JACC = os.path.join(ROOT, "results/four_catalog_jaccard.tsv")

ARMS = ["herptile", "ehi", "youngblut", "gtdb_ref"]
N_RESAMPLE = 499
random.seed(20260804)


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


# ------------------------------------------------------------ arm + genus
arm_of, gid_of = {}, {}
with open(MAP) as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        arm_of[f[0]] = f[1]
        gid_of[f[0]] = f[2]

genus = {}
for r in read_tsv(HERP_SUM):
    genus[("herptile", r["user_genome"].strip())] = norm(
        parse_tax(r["classification"]).get("g", ""))
for r in read_tsv(EHI_SUM):
    genus[("ehi", r["user_genome"].strip())] = norm(
        parse_tax(r["classification"]).get("g", ""))
for r in read_tsv(YB_SUM):
    genus[("youngblut", r["user_genome"].strip())] = norm(
        parse_tax(r["classification"]).get("g", ""))

n_ref = 0
with open(GTDB_TAX) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        acc = f[0].strip().replace("GB_", "").replace("RS_", "")
        t = parse_tax(f[1])
        if t.get("f", "") == "Ruminococcaceae":
            genus[("gtdb_ref", acc)] = norm(t.get("g", ""))
            n_ref += 1
print("genus labels loaded: %d (gtdb_ref %d)" % (len(genus), n_ref))

# ------------------------------------------------------------ clusters
clusters = defaultdict(list)
with open(CDB) as fh:
    hdr = fh.readline().rstrip("\n").replace('"', "").split(",")
    gi, si = hdr.index("genome"), hdr.index("secondary_cluster")
    for line in fh:
        f = line.rstrip("\n").replace('"', "").split(",")
        if len(f) > max(gi, si):
            clusters[f[si]].append(f[gi])
print("clusters: %d" % len(clusters))

# one unit per cluster per arm: an SGB counts once for each arm present
unit_genus = defaultdict(list)   # arm -> list of genus labels, one per cluster
unassigned = Counter()
for cid, gs in clusters.items():
    per_arm = defaultdict(list)
    for g in gs:
        a = arm_of.get(g, "?")
        lab = genus.get((a, gid_of.get(g, "")), None)
        per_arm[a].append(lab)
    for a, labs in per_arm.items():
        named = [x for x in labs if x]
        if named:
            unit_genus[a].append(Counter(named).most_common(1)[0][0])
        else:
            unassigned[a] += 1

print()
print("=" * 74)
print("UNITS: ONE CROSS-CATALOG SGB PER ARM PER CLUSTER")
print("=" * 74)
print("  %-10s %8s %10s %10s" % ("arm", "units", "named", "no genus"))
sets = {}
for a in ARMS:
    labs = unit_genus.get(a, [])
    sets[a] = set(labs)
    print("  %-10s %8d %10d %10d"
          % (a, len(labs) + unassigned.get(a, 0), len(labs), unassigned.get(a, 0)))
print("  Every arm is now the same unit (95%% ANI cluster), so counts ARE")
print("  comparable across columns. This was not true before dereplication.")

# ------------------------------------------------------------ table
counts = {a: Counter(unit_genus.get(a, [])) for a in ARMS}
allg = set()
for a in ARMS:
    allg |= sets[a]

ENDO = ("ehi", "youngblut")


def block_of(g):
    h = g in sets["herptile"]
    e = any(g in sets[a] for a in ENDO)
    r = g in sets["gtdb_ref"]
    if h and not e and not r:
        return "1_herptile_only"
    if h and not e and r:
        return "2_herptile_and_reference_only"
    if h and e:
        return "3_herptile_and_endotherm"
    return "4_no_herptile"


# V2 LABELS. block_of() tests ONLY EHI mammals and Youngblut, so the label
# must name those two catalogs and must not generalise to endotherm guts.
BLOCK_LABEL = {
    "1_herptile_only": "herptile only, not recovered by EHI, Youngblut or GTDB references",
    "2_herptile_and_reference_only": "herptile plus GTDB reference, not recovered by EHI or Youngblut",
    "3_herptile_and_endotherm": "herptile and at least one of EHI or Youngblut",
    "4_no_herptile": "absent from herptile",
}

blocks = defaultdict(list)
for g in allg:
    blocks[block_of(g)].append(g)

print()
print("=" * 74)
print("BLOCKS")
print("=" * 74)
for b in sorted(BLOCK_LABEL):
    gl = blocks.get(b, [])
    units = sum(counts["herptile"].get(g, 0) for g in gl)
    print("  %-34s %3d genera" % (b, len(gl)))
    print("      %s" % BLOCK_LABEL[b])
    if b != "4_no_herptile":
        print("      holding %d herptile SGBs" % units)

print()
print("  WORDING, MANDATORY: block 2 is NOT 'absent from endotherm guts'.")
print("  Those genera hold 62 GTDB reference genomes and 42 are gut or")
print("  faecal, including 17 ruminant GI tract. RGIG3102 alone is 20 of the")
print("  62 and nearly all ruminant. This is DIFFERENTIAL RECOVERY by two")
print("  named catalogs, not compositional absence.")
print("  file: results/amphibian_only_genus_sources.tsv")
print()
print("  ARMS: 'ehi' is EHI MAMMAL/BIRD only. EHI newts are not tested here.")
print("  For the newt comparison see results/amphibian_genus_replication.tsv")
print("  (in_both_amphibian x in_any_endotherm = 8 / 6 / 25 / 98).")

with open(OUT_TABLE, "w") as f:
    f.write("genus\tblock\tblock_label\therptile\tehi\tyoungblut\tgtdb_ref\n")
    for b in sorted(BLOCK_LABEL):
        for g in sorted(blocks.get(b, []),
                        key=lambda x: -(counts["herptile"].get(x, 0)
                                        + counts["ehi"].get(x, 0)
                                        + counts["youngblut"].get(x, 0))):
            f.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\n"
                    % (g, b, BLOCK_LABEL[b],
                       counts["herptile"].get(g, 0), counts["ehi"].get(g, 0),
                       counts["youngblut"].get(g, 0), counts["gtdb_ref"].get(g, 0)))
print()
print("  wrote %s" % OUT_TABLE)

# ------------------------------------------------------------ jaccard
def jac(a, b):
    u = a | b
    return len(a & b) / float(len(u)) if u else 0.0


def q(s, p):
    return s[int(round(p * (len(s) - 1)))] if s else float("nan")


print()
print("=" * 74)
print("JACCARD ON DEREPLICATED UNITS")
print("=" * 74)
PAIRS = [("herptile", "ehi"), ("herptile", "youngblut"), ("herptile", "gtdb_ref"),
         ("youngblut", "ehi"), ("gtdb_ref", "ehi"), ("gtdb_ref", "youngblut")]
print("  %-26s %8s %8s %8s %7s" % ("pair", "A gen", "B gen", "shared", "J"))
for a, b in PAIRS:
    print("  %-26s %8d %8d %8d %7.3f"
          % (a + " vs " + b, len(sets[a]), len(sets[b]),
             len(sets[a] & sets[b]), jac(sets[a], sets[b])))

kmin = min(len(sets[a]) for a in ARMS)
print()
print("  RICHNESS-MATCHED to k = %d (%d resamples)" % (kmin, N_RESAMPLE))
out = open(OUT_JACC, "w")
out.write("pair\tn_a\tn_b\tshared\tjaccard_raw\tk\tjaccard_matched\tlo\thi\n")
for a, b in PAIRS:
    la, lb = list(sets[a]), list(sets[b])
    js = [jac(set(random.sample(la, kmin)), set(random.sample(lb, kmin)))
          for _ in range(N_RESAMPLE)]
    s = sorted(js)
    m = sum(js) / float(len(js))
    print("  %-26s %.3f [%.3f, %.3f]" % (a + " vs " + b, m, q(s, 0.025), q(s, 0.975)))
    out.write("%s_vs_%s\t%d\t%d\t%d\t%.4f\t%d\t%.4f\t%.4f\t%.4f\n"
              % (a, b, len(sets[a]), len(sets[b]), len(sets[a] & sets[b]),
                 jac(sets[a], sets[b]), kmin, m, q(s, 0.025), q(s, 0.975)))
out.close()
print()
print("  wrote %s" % OUT_JACC)
print()
print("  Compare herptile-vs-ehi against youngblut-vs-ehi. Before")
print("  dereplication that contrast was 0.084 vs 0.314. If it holds now that")
print("  ehi is 80 units rather than 280 genomes, redundancy was not driving it.")
print("FOUR_CATALOG_GENUS_TABLE_V2_20260806_COMPLETE")
