# Reference antiSMASH output is parsed. The manifest schema differs from the amphibian one.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/parse_antismash_refs.py
# Output: results/bgc_regions_refs.tsv, results/bgc_per_genome_refs.tsv
# Reference-arm antiSMASH parser. Separate from scripts/parse_antismash.py,
# which is hard-coded to work/bgc and results/bgc_*.tsv and would overwrite
# the amphibian R9 tables.
#
# Manifest schema differs from the amphibian one: columns are
# index, safe_name, accession, genus, staged_path, source_path.
# There is no arm, family or host column; all 1,247 are GTDB r220
# Ruminococcaceae references, and ncbi_host_name is blank for these
# genomes anyway (see R7).
#
# Outputs are new files. Nothing existing is touched.

import os, re
from collections import Counter, defaultdict

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
WORK = os.path.join(ROOT, "work/bgc_refs")
ASDIR = os.path.join(WORK, "antismash")
MAN = os.path.join(WORK, "bgc_ref_manifest.tsv")
OUT_REG = os.path.join(ROOT, "results/bgc_regions_refs.tsv")
OUT_GEN = os.path.join(ROOT, "results/bgc_per_genome_refs.tsv")

for p in (OUT_REG, OUT_GEN):
    if os.path.exists(p):
        raise SystemExit("REFUSING TO OVERWRITE %s, move it first" % p)

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
print("  distinct genera  : %d" % len(set(r["genus"] for r in man.values())))

if not os.path.isdir(ASDIR):
    raise SystemExit("no antismash output at %s" % ASDIR)
dirs = set(e.name for e in os.scandir(ASDIR) if e.is_dir())
print("  output directories: %d" % len(dirs))

# A directory is only evidence of a finished run if it holds the json.
# One genome (gtdbref__GCA_000174895.1) previously had a directory from an
# aborted duplicate task; counting directories would have scored it zero.
have_json = set()
for d in dirs:
    if os.path.exists(os.path.join(ASDIR, d, d + ".json")):
        have_json.add(d)
print("  directories with a json: %d" % len(have_json))

missing = sorted(set(man) - have_json)
print("  manifest genomes with NO json: %d" % len(missing))
for s in missing[:10]:
    print("     %s" % s)
extra = sorted(have_json - set(man))
if extra:
    print("  json present but NOT in manifest: %d" % len(extra))
    for s in extra[:10]:
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
            genome=safe,
            region=g,
            product=prods[0] if prods else "unknown",
            edge=(edges[0].lower() == "true") if edges else None,
        ))

n_edge = sum(1 for r in regions if r["edge"] is True)
n_comp = sum(1 for r in regions if r["edge"] is False)
n_unk = sum(1 for r in regions if r["edge"] is None)

print()
print("=" * 74)
print("REGIONS")
print("=" * 74)
print("  genomes parsed        : %d" % len(have_json & set(man)))
print("  total regions         : %d" % len(regions))
print("  complete              : %d" % n_comp)
print("  on a contig edge      : %d (%.1f%%)"
      % (n_edge, 100.0 * n_edge / len(regions) if regions else 0))
print("  contig_edge not parsed: %d" % n_unk)

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
vals_all = [per[g]["all"] for g in parsed]
vals_c = [per[g]["comp"] for g in parsed]
n = float(len(parsed)) if parsed else 1.0

print()
print("=" * 74)
print("PER GENOME, WHOLE REFERENCE ARM")
print("=" * 74)
print("  genomes            : %d" % len(parsed))
print("  mean regions (all) : %.2f" % (sum(vals_all) / n))
print("  mean COMPLETE      : %.2f" % (sum(vals_c) / n))
print("  %% with any region  : %.1f%%" % (100.0 * sum(1 for v in vals_all if v > 0) / n))
print("  %% with a complete  : %.1f%%" % (100.0 * sum(1 for v in vals_c if v > 0) / n))
print()
print("  The amphibian arm gave mean COMPLETE 0.48 (EHI newt), 0.50")
print("  (herptile amphibian), 0.61 (herptile reptile). Compare to the")
print("  mean COMPLETE above, never to the all-region totals: 91.9% of GTDB")
print("  r220 Ruminococcaceae are themselves MAGs (R7) and fragmentation")
print("  depresses both arms, but not necessarily equally.")

print()
print("=" * 74)
print("BY GENUS, GENERA WITH >= 10 GENOMES")
print("=" * 74)
bygen = defaultdict(list)
for g in parsed:
    bygen[man[g]["genus"]].append(g)
print("  %-28s %6s %10s %12s" % ("genus", "n", "mean all", "mean complete"))
for gen in sorted(bygen, key=lambda k: -len(bygen[k])):
    gs = bygen[gen]
    if len(gs) < 10:
        continue
    print("  %-28s %6d %10.2f %12.2f"
          % (gen[:28], len(gs),
             sum(per[g]["all"] for g in gs) / float(len(gs)),
             sum(per[g]["comp"] for g in gs) / float(len(gs))))

with open(OUT_REG, "w") as f:
    f.write("genome\tarm\taccession\tgenus\tregion\tproduct\tcontig_edge\n")
    for r in regions:
        m = man.get(r["genome"], {})
        f.write("%s\treference\t%s\t%s\t%s\t%s\t%s\n"
                % (r["genome"], m.get("accession", ""), m.get("genus", ""),
                   r["region"], r["product"],
                   "" if r["edge"] is None else ("yes" if r["edge"] else "no")))

with open(OUT_GEN, "w") as f:
    f.write("genome\tarm\taccession\tgenus\tn_regions\tn_complete\tn_edge\tproducts\n")
    for safe in parsed:
        m = man[safe]
        d = per.get(safe) or dict(all=0, comp=0, edge=0, prods=Counter())
        f.write("%s\treference\t%s\t%s\t%d\t%d\t%d\t%s\n"
                % (safe, m["accession"], m["genus"],
                   d["all"], d["comp"], d["edge"],
                   ";".join("%s:%d" % kv for kv in sorted(d["prods"].items()))))

print()
print("wrote %s" % OUT_REG)
print("wrote %s" % OUT_GEN)
print()
print("  Genomes with no json are EXCLUDED from both files and from every")
print("  mean above. The denominator is the parsed count, not 1,247.")
print("PARSE_ANTISMASH_REFS_V1_20260805_COMPLETE")
