# Amphibian antiSMASH output is parsed into region and per-genome tables.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/parse_antismash.py
# Output: results/bgc_regions.tsv, results/bgc_per_genome.tsv
# PARSE_ANTISMASH_V1_20260805
import os, re, json, gzip
from collections import Counter, defaultdict

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
WORK = os.path.join(ROOT, "work/bgc")
ASDIR = os.path.join(WORK, "antismash")
MAN = os.path.join(WORK, "bgc_manifest.tsv")
OUT_REG = os.path.join(ROOT, "results/bgc_regions.tsv")
OUT_GEN = os.path.join(ROOT, "results/bgc_per_genome.tsv")


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
print("  manifest entries: %d" % len(man))
print("  by arm : %s" % dict(Counter(r["arm"] for r in man.values())))
print("  by host: %s" % dict(Counter(r["host"] for r in man.values())))

if not os.path.isdir(ASDIR):
    raise SystemExit("no antismash output at %s" % ASDIR)
dirs = [e.name for e in os.scandir(ASDIR) if e.is_dir()]
print("  antismash output directories: %d" % len(dirs))

# ------------------------------------------------------ inspect layout
sample = sorted(dirs)[0] if dirs else None
if sample:
    files = sorted(os.listdir(os.path.join(ASDIR, sample)))
    print()
    print("  files in %s:" % sample)
    for f in files[:14]:
        print("     %s" % f)
    if len(files) > 14:
        print("     ... and %d more" % (len(files) - 14))

# ------------------------------------------------------ parse regions
# Region GBK files carry /product= and /contig_edge= in the region feature.
# contig_edge=True means the cluster runs off the end of a contig and is
# therefore incomplete. Report complete and edge counts SEPARATELY: MAG BGC
# totals are not comparable to isolate genomes because of fragmentation.
prod_re = re.compile(r'/product="([^"]+)"')
edge_re = re.compile(r'/contig_edge="([^"]+)"')

regions = []
no_output = []
for safe in sorted(man):
    d = os.path.join(ASDIR, safe)
    if not os.path.isdir(d):
        no_output.append(safe)
        continue
    gbks = [f for f in os.listdir(d) if ".region" in f and f.endswith(".gbk")]
    for g in sorted(gbks):
        p = os.path.join(d, g)
        try:
            with open(p) as fh:
                head = fh.read(20000)
        except Exception:
            continue
        prods = prod_re.findall(head)
        edges = edge_re.findall(head)
        prod = prods[0] if prods else "unknown"
        edge = (edges[0].lower() == "true") if edges else None
        regions.append(dict(genome=safe, region=g, product=prod, edge=edge))

print()
print("=" * 74)
print("REGIONS")
print("=" * 74)
print("  genomes with no output directory: %d" % len(no_output))
for s in no_output[:5]:
    print("     %s" % s)
print("  total regions parsed: %d" % len(regions))
n_edge = sum(1 for r in regions if r["edge"] is True)
n_comp = sum(1 for r in regions if r["edge"] is False)
n_unk = sum(1 for r in regions if r["edge"] is None)
print("  complete            : %d" % n_comp)
print("  on a contig edge    : %d (%.1f%%)"
      % (n_edge, 100.0 * n_edge / len(regions) if regions else 0))
print("  contig_edge not parsed: %d" % n_unk)
print()
print("  A cluster on a contig edge is truncated by assembly, not biology.")
print("  BiG-SLiCE builds models from non-fragmented BGCs only and maps")
print("  partials back afterwards, which is the workflow to follow.")

# ------------------------------------------------------ by product
print()
print("=" * 74)
print("PRODUCT CLASSES")
print("=" * 74)
pc_all = Counter(r["product"] for r in regions)
pc_comp = Counter(r["product"] for r in regions if r["edge"] is False)
print("  %-34s %10s %10s" % ("product", "all", "complete"))
for p, n in pc_all.most_common(20):
    print("  %-34s %10d %10d" % (p[:34], n, pc_comp.get(p, 0)))
if len(pc_all) > 20:
    print("  ... and %d more classes" % (len(pc_all) - 20))

# ------------------------------------------------------ per genome
per = defaultdict(lambda: dict(all=0, comp=0, edge=0, prods=Counter()))
for r in regions:
    d = per[r["genome"]]
    d["all"] += 1
    if r["edge"] is True:
        d["edge"] += 1
    elif r["edge"] is False:
        d["comp"] += 1
    d["prods"][r["product"]] += 1

print()
print("=" * 74)
print("PER GENOME, BY ARM AND HOST")
print("=" * 74)
groups = defaultdict(list)
for safe, r in man.items():
    if safe in per or safe not in no_output:
        groups[(r["arm"], r["host"])].append(safe)

print("  %-16s %-12s %6s %10s %12s %12s"
      % ("arm", "host", "n", "mean all", "mean complete", "% with any"))
for k in sorted(groups):
    gs = [g for g in groups[k] if g in per or g not in no_output]
    vals_all = [per[g]["all"] for g in gs]
    vals_c = [per[g]["comp"] for g in gs]
    withany = sum(1 for v in vals_all if v > 0)
    if not gs:
        continue
    print("  %-16s %-12s %6d %10.2f %12.2f %11.1f%%"
          % (k[0], k[1], len(gs),
             sum(vals_all) / float(len(gs)), sum(vals_c) / float(len(gs)),
             100.0 * withany / len(gs)))

# ------------------------------------------------------ MIBiG
print()
print("=" * 74)
print("MIBiG COMPARISON (knownclusterblast)")
print("=" * 74)
kcb_dirs = 0
kcb_hits = 0
for safe in sorted(per):
    kd = os.path.join(ASDIR, safe, "knownclusterblast")
    if os.path.isdir(kd):
        kcb_dirs += 1
        for f in os.listdir(kd):
            if f.endswith(".txt"):
                try:
                    with open(os.path.join(kd, f)) as fh:
                        t = fh.read()
                    if "Significant hits" in t and "BGC" in t:
                        kcb_hits += 1
                except Exception:
                    pass
print("  genomes with a knownclusterblast directory: %d" % kcb_dirs)
print("  region files reporting a significant MIBiG hit: %d" % kcb_hits)
if kcb_dirs == 0:
    print("  none found; --cb-knownclusters may not have produced output.")
print()
print("  For context, Youngblut et al. 2020 mSystems reported 1,986 BGCs of")
print("  which only 23 clustered with any MIBiG reference. A low match rate")
print("  here is expected and is itself the novelty statement.")

# ------------------------------------------------------ write
with open(OUT_REG, "w") as f:
    f.write("genome\tarm\thost\tfamily\tgenus\tregion\tproduct\tcontig_edge\n")
    for r in regions:
        m = man.get(r["genome"], {})
        f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                % (r["genome"], m.get("arm", ""), m.get("host", ""),
                   m.get("family", ""), m.get("genus", ""),
                   r["region"], r["product"],
                   "" if r["edge"] is None else ("yes" if r["edge"] else "no")))
with open(OUT_GEN, "w") as f:
    f.write("genome\tarm\thost\tfamily\tgenus\tn_regions\tn_complete\tn_edge\tproducts\n")
    for safe in sorted(man):
        m = man[safe]
        d = per.get(safe)
        if d is None and safe in no_output:
            continue
        d = d or dict(all=0, comp=0, edge=0, prods=Counter())
        f.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\n"
                % (safe, m["arm"], m["host"], m["family"], m["genus"],
                   d["all"], d["comp"], d["edge"],
                   ";".join("%s:%d" % kv for kv in sorted(d["prods"].items()))))
print()
print("wrote %s" % OUT_REG)
print("wrote %s" % OUT_GEN)

print()
print("=" * 74)
print("WHAT IS AND IS NOT CLAIMABLE")
print("=" * 74)
print("  Claimable: the number and classes of BGC regions detected in these")
print("  genomes, reported separately for complete and contig-edge regions.")
print("  NOT claimable without more work: comparison of these totals to")
print("  isolate genomes (fragmentation depresses MAG counts), novelty")
print("  (needs BiG-SCAPE or BiG-SLiCE clustering into gene cluster families")
print("  and comparison to MIBiG AND BiG-FAM), or any link to host biology.")

# SENTINEL_END
