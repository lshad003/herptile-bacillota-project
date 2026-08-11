# Endotherm antiSMASH output is parsed, with the two source catalogs reported separately.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/parse_antismash_endo.py
# Output: results/bgc_regions_endo.tsv, results/bgc_per_genome_endo.tsv
# PARSE_ANTISMASH_ENDO_V1_20260806
# Endotherm-arm antiSMASH parser, job 27253466.
# 280 EHI mammal + 48 Youngblut Ruminococcaceae, reported SEPARATELY.
#
# Separate from parse_antismash.py (hard-coded to work/bgc) and
# parse_antismash_refs.py (work/bgc_refs). Writes NEW files; nothing existing
# is touched.
#
# WHY THE EHI MAMMAL ARM IS THE POINT: EHI mammals and EHI newts come from
# one consortium, one pipeline, one data release, so assembly contiguity is
# matched BY CONSTRUCTION rather than by luck. R9 established that
# complete-BGC recovery is a function of N50 that does not differ by catalog
# (at matched N50: 0.34/0.14/0.37 in the 20-40 kb bin, 0.67/0.65/0.68 in
# 40-80 kb). This arm tests whether it also does not differ by HOST.
#
# NAME THE ARM BY WHAT IT CONTAINS. EHI mammal is Rodentia 119, Carnivora 82,
# Lagomorpha 76, Diprotodontia 2, Psittaciformes 1. NO ruminants, NO humans.
# It is not "endotherm gut".
#
# YOUNGBLUT IS SECONDARY: 79.2% already sit inside a GTDB r220 species
# cluster, so the arm is not independent of the reference set, and 22 of 48
# are chicken. Do not build a claim on it alone.
#
# A DIRECTORY IS NOT EVIDENCE OF A FINISHED RUN. Count the json.

import os, re, sys
from collections import Counter, defaultdict

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
WORK = os.path.join(ROOT, "work/bgc_endo")
ASDIR = os.path.join(WORK, "antismash")
MAN = os.path.join(WORK, "bgc_endo_manifest.tsv")
OUT_REG = os.path.join(ROOT, "results/bgc_regions_endo.tsv")
OUT_GEN = os.path.join(ROOT, "results/bgc_per_genome_endo.tsv")

for p in (OUT_REG, OUT_GEN):
    if os.path.exists(p):
        raise SystemExit("REFUSING TO OVERWRITE %s, move it first" % p)
for p in (MAN,):
    if not os.path.exists(p):
        sys.exit("MISSING: %s" % p)
if not os.path.isdir(ASDIR):
    sys.exit("no antismash output at %s" % ASDIR)


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


man = {r["safe_name"]: r for r in read_tsv(MAN)}
print("=" * 74)
print("INPUT")
print("=" * 74)
print("  manifest entries : %d" % len(man))
print("  by arm           : %s" % dict(Counter(r["arm"] for r in man.values())))

dirs = set(e.name for e in os.scandir(ASDIR) if e.is_dir())
have_json = set(d for d in dirs
                if os.path.exists(os.path.join(ASDIR, d, d + ".json")))
print("  output dirs      : %d" % len(dirs))
print("  dirs with a json : %d" % len(have_json))
missing = sorted(set(man) - have_json)
print("  manifest genomes with NO json: %d" % len(missing))
for s in missing[:10]:
    print("     %s" % s)

prod_re = re.compile(r'/product="([^"]+)"')
edge_re = re.compile(r'/contig_edge="([^"]+)"')

regions = []
for safe in sorted(man):
    if safe not in have_json:
        continue
    d = os.path.join(ASDIR, safe)
    gbks = [f for f in os.listdir(d) if ".region" in f and f.endswith(".gbk")]
    for g in sorted(gbks):
        try:
            with open(os.path.join(d, g)) as fh:
                head = fh.read(20000)
        except Exception:
            continue
        prods = prod_re.findall(head)
        edges = edge_re.findall(head)
        regions.append(dict(
            genome=safe, arm=man[safe]["arm"], region=g,
            product=prods[0] if prods else "unknown",
            edge=(edges[0].lower() == "true") if edges else None))

per = defaultdict(lambda: dict(all=0, comp=0, edge=0, prods=Counter()))
for r in regions:
    d = per[r["genome"]]
    d["all"] += 1
    if r["edge"] is True:
        d["edge"] += 1
    elif r["edge"] is False:
        d["comp"] += 1
    d["prods"][r["product"]] += 1

parsed = sorted(have_json & set(man))
ARMS = ["ehi_mammal", "youngblut"]

print()
print("=" * 74)
print("REGIONS, BY ARM")
print("=" * 74)
print("  %-14s %6s %9s %9s %9s %8s" %
      ("arm", "n", "regions", "complete", "edge", "edge %"))
for a in ARMS:
    gs = [g for g in parsed if man[g]["arm"] == a]
    if not gs:
        continue
    tot = sum(per[g]["all"] for g in gs)
    c = sum(per[g]["comp"] for g in gs)
    e = sum(per[g]["edge"] for g in gs)
    print("  %-14s %6d %9d %9d %9d %7.1f%%"
          % (a, len(gs), tot, c, e, 100.0 * e / tot if tot else 0.0))
unk = sum(1 for r in regions if r["edge"] is None)
print("  contig_edge not parsed: %d" % unk)

print()
print("=" * 74)
print("PER GENOME, BY ARM")
print("=" * 74)
print("  %-14s %6s %12s %14s %12s %12s" %
      ("arm", "n", "mean all", "mean COMPLETE", "% any", "% any compl"))
for a in ARMS:
    gs = [g for g in parsed if man[g]["arm"] == a]
    if not gs:
        continue
    n = float(len(gs))
    va = [per[g]["all"] for g in gs]
    vc = [per[g]["comp"] for g in gs]
    print("  %-14s %6d %12.2f %14.2f %11.1f%% %11.1f%%"
          % (a, len(gs), sum(va) / n, sum(vc) / n,
             100.0 * sum(1 for v in va if v > 0) / n,
             100.0 * sum(1 for v in vc if v > 0) / n))
print()
print("  For comparison, Ruminococcaceae-only mean COMPLETE from R9:")
print("    UHM amphibian 0.28 (n=216) | EHI newt 0.33 (n=162)")
print("    GTDB reference 0.64 (n=1247)")
print("  THE COMPARISON THAT MATTERS is EHI mammal vs EHI newt: same")
print("  consortium, pipeline and data release, so contiguity is matched by")
print("  construction. Do NOT compare either to GTDB, which fails every")
print("  assembly-quality balance test (N50 SMD 0.31, contig count 0.76).")
print("  Assembly quality for this arm is NOT yet measured; run it before")
print("  any comparison and add these 328 genomes to")
print("  results/assembly_quality_arms.tsv.")

print()
print("=" * 74)
print("EHI MAMMAL, BY HOST ORDER")
print("=" * 74)
byho = defaultdict(list)
for g in parsed:
    if man[g]["arm"] == "ehi_mammal":
        byho[man[g].get("host_order", "") or "(blank)"].append(g)
print("  %-18s %6s %14s %12s" % ("host order", "n", "mean COMPLETE", "mean all"))
for ho in sorted(byho, key=lambda k: -len(byho[k])):
    gs = byho[ho]
    n = float(len(gs))
    print("  %-18s %6d %14.2f %12.2f"
          % (ho, len(gs), sum(per[g]["comp"] for g in gs) / n,
             sum(per[g]["all"] for g in gs) / n))
print("  Orders with n < 10 are shown for completeness only.")

print()
print("=" * 74)
print("PRODUCT CLASSES, COMPLETE REGIONS ONLY")
print("=" * 74)
pc = defaultdict(Counter)
for r in regions:
    if r["edge"] is False:
        pc[r["arm"]][r["product"]] += 1
allp = Counter()
for a in ARMS:
    allp.update(pc[a])
print("  %-34s %12s %12s" % ("product", "ehi_mammal", "youngblut"))
for p, _ in allp.most_common(12):
    print("  %-34s %12d %12d"
          % (p[:34], pc["ehi_mammal"].get(p, 0), pc["youngblut"].get(p, 0)))

with open(OUT_REG, "w") as f:
    f.write("genome\tarm\tgenome_id\tgenus\thost_order\tregion\tproduct\tcontig_edge\n")
    for r in regions:
        m = man.get(r["genome"], {})
        f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                % (r["genome"], r["arm"], m.get("genome_id", ""),
                   m.get("genus", ""), m.get("host_order", ""),
                   r["region"], r["product"],
                   "" if r["edge"] is None else ("yes" if r["edge"] else "no")))

with open(OUT_GEN, "w") as f:
    f.write("genome\tarm\tgenome_id\tgenus\thost_species\thost_order\t"
            "completeness\tcontigs\tn_regions\tn_complete\tn_edge\tproducts\n")
    for safe in parsed:
        m = man[safe]
        d = per.get(safe) or dict(all=0, comp=0, edge=0, prods=Counter())
        f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\n"
                % (safe, m["arm"], m["genome_id"], m["genus"],
                   m.get("host_species", ""), m.get("host_order", ""),
                   m.get("completeness", ""), m.get("contigs", ""),
                   d["all"], d["comp"], d["edge"],
                   ";".join("%s:%d" % kv for kv in sorted(d["prods"].items()))))

print()
print("wrote %s" % OUT_REG)
print("wrote %s" % OUT_GEN)
print()
print("CAPTION AND WRITING REQUIREMENTS:")
print("  Ruminococcaceae only. Complete and contig-edge separate throughout.")
print("  Call the arm 'EHI mammal', never 'endotherm gut': it is 277 of 280")
print("  Rodentia, Carnivora and Lagomorpha, with no ruminants and no humans.")
print("  antiSMASH 7.1.0, --genefinding-tool prodigal-m, --cb-knownclusters.")
print("  Exploratory, not preregistered. No test has been run.")
print("PARSE_ANTISMASH_ENDO_V1_20260806_COMPLETE")
