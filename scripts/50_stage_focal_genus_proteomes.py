# Proteomes are staged for the two interleaved genera.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/stage_focal_genus_proteomes.py
# Output: work/focal_genus_pangenome/proteomes/, focal_genome_manifest.tsv
# STAGE_FOCAL_GENUS_V1_20260804
import os, gzip, shutil
from collections import Counter, defaultdict

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
META = os.path.join(ROOT, "work/rep_tree/figure_tree_metadata.tsv")
GTDB_TAX = "/srv/projects/db/gtdbtk/220/taxonomy/gtdb_taxonomy.tsv"
REFPATHS = os.path.join(ROOT, "data/tasks/gtdb_ruminococcaceae_paths.tsv")
DREP_MAP = os.path.join(ROOT, "work/pooled_drep_rum/genome_arms.tsv")
DREP_GEN = os.path.join(ROOT, "work/pooled_drep_rum/genomes")

WORK = os.path.join(ROOT, "work/focal_genus_pangenome")
GDIR = os.path.join(WORK, "genomes")
OUT_MAN = os.path.join(WORK, "focal_genome_manifest.tsv")

FOCAL = ["Anaerotruncus", "UBA866"]
AMPH_ARMS = {"herptile", "ehi_amphibian"}


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


for d in (WORK, GDIR):
    if not os.path.isdir(d):
        os.makedirs(d)

# reference genus map, accession is the LAST pipe field of the tip label
ref_genus = {}
with open(GTDB_TAX) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        t = parse_tax(f[1])
        if t.get("f", "") == "Ruminococcaceae":
            ref_genus[f[0].strip().replace("GB_", "").replace("RS_", "")] = t.get("g", "")

meta = read_tsv(META)
for m in meta:
    if m["arm"] == "gtdb_ref" and not m["genus"]:
        m["genus"] = ref_genus.get(m["tip"].split("|")[-1], "")

targets = [m for m in meta if m["genus"] in FOCAL and m["arm"] != "outgroup"]
print("=" * 74)
print("FOCAL GENERA: %s" % ", ".join(FOCAL))
print("=" * 74)
byg = defaultdict(Counter)
for m in targets:
    grp = "amphibian" if m["arm"] in AMPH_ARMS else (
        "reference" if m["arm"] == "gtdb_ref" else "endotherm")
    byg[m["genus"]][grp] += 1
for g in FOCAL:
    d = byg[g]
    print("  %-16s amphibian %3d | endotherm %3d | reference %3d | total %3d"
          % (g, d["amphibian"], d["endotherm"], d["reference"], sum(d.values())))
n_amph = sum(1 for m in targets if m["arm"] in AMPH_ARMS)
print("  POOLED: %d amphibian vs %d non-amphibian, %d genomes total"
      % (n_amph, len(targets) - n_amph, len(targets)))

# ------------------------------------------------------- resolve fastas
drep_path = {}
with open(DREP_MAP) as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        drep_path[(f[1], f[2])] = os.path.join(DREP_GEN, f[0])

ref_path = {}
for r in read_tsv(REFPATHS):
    ref_path[r["accession"].strip()] = r["path"].strip()

print()
print("=" * 74)
print("RESOLVING GENOME FILES")
print("=" * 74)
resolved, missing = [], []
for m in targets:
    arm, gid = m["arm"], m["genome"]
    if arm == "gtdb_ref":
        acc = m["tip"].split("|")[-1]
        p = ref_path.get(acc, "")
        gid = acc
    else:
        p = drep_path.get((arm, gid), "")
    if p and os.path.exists(p) and os.path.getsize(p) > 0:
        resolved.append((m, gid, p))
    else:
        missing.append((arm, gid, p))
print("  resolved: %d, missing: %d" % (len(resolved), len(missing)))
for a, g, p in missing[:6]:
    print("     %-14s %-24s %s" % (a, g, p if p else "<no path>"))

# ------------------------------------------------------------ stage
print()
print("=" * 74)
print("STAGING UNCOMPRESSED NUCLEOTIDE FASTAS")
print("=" * 74)
n_new = n_skip = n_fail = 0
staged = []
for m, gid, src in resolved:
    grp = "amphibian" if m["arm"] in AMPH_ARMS else (
        "reference" if m["arm"] == "gtdb_ref" else "endotherm")
    safe = "%s__%s__%s" % (m["genus"], grp, gid)
    safe = safe.replace("/", "_").replace(" ", "_")
    dst = os.path.join(GDIR, safe + ".fna")
    if os.path.exists(dst) and os.path.getsize(dst) > 0:
        staged.append((m, gid, safe, grp))
        n_skip += 1
        continue
    try:
        if src.endswith(".gz"):
            with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
                shutil.copyfileobj(fi, fo)
        else:
            shutil.copyfile(src, dst)
        if os.path.getsize(dst) == 0:
            raise IOError("zero bytes")
        staged.append((m, gid, safe, grp))
        n_new += 1
    except Exception as ex:
        n_fail += 1
        if n_fail <= 5:
            print("     FAIL %s: %s" % (gid, ex))
        if os.path.exists(dst):
            os.remove(dst)
print("  staged %d, already present %d, failed %d" % (n_new, n_skip, n_fail))

files = [e for e in os.scandir(GDIR) if e.name.endswith(".fna")]
tot = sum(e.stat().st_size for e in files)
print("  files: %d, %.1f MB" % (len(files), tot / 1e6))
print("  zero-byte: %d" % sum(1 for e in files if e.stat().st_size == 0))

with open(OUT_MAN, "w") as f:
    f.write("staged_name\tgenus\tgroup\tarm\tgenome\thost_detail\tgap_fraction\n")
    for m, gid, safe, grp in sorted(staged, key=lambda x: x[2]):
        f.write("%s.fna\t%s\t%s\t%s\t%s\t%s\t%s\n"
                % (safe, m["genus"], grp, m["arm"], gid,
                   m["host_detail"], m["gap_fraction"]))
print("  wrote %s (%d rows)" % (OUT_MAN, len(staged)))

print()
print("  final composition:")
c = Counter((m["genus"], grp) for m, _, _, grp in staged)
for k in sorted(c):
    print("     %-16s %-12s %3d" % (k[0], k[1], c[k]))
print()
print("  Design: %d amphibian vs %d non-amphibian, pooled across two genera,"
      % (sum(1 for m, _, _, g in staged if g == "amphibian"),
         sum(1 for m, _, _, g in staged if g != "amphibian")))
print("  with genus as a covariate. Chosen because these are the only two")
print("  genera where amphibian and non-amphibian genomes interleave on the")
print("  tree (Anaerotruncus 7 transitions, UBA866 5), so a host effect is")
print("  not simply a single lineage effect.")

# SENTINEL_END
