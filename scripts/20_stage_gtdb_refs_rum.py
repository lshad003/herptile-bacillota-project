# GTDB r220 Ruminococcaceae reference genomes are added to the dereplication input set.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/stage_gtdb_refs_rum.py
# Output: work/pooled_drep_rum/genomes/, genome_arms.tsv, genome_info.csv
# STAGE_GTDB_REFS_RUM_V1_20260804
import os, gzip, shutil
from collections import Counter

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
WORK = os.path.join(ROOT, "work/pooled_drep_rum")
GDIR = os.path.join(WORK, "genomes")
MAP = os.path.join(WORK, "genome_arms.tsv")
INFO = os.path.join(WORK, "genome_info.csv")
REFS = os.path.join(ROOT, "data/tasks/gtdb_ruminococcaceae_paths.tsv")
META = "/srv/projects/db/gtdbtk/220/metadata/gtdb_release_metadata.tsv"

# GTDB genomes carry their own completeness/contamination in the release
# metadata. If that file is absent we fall back to a conservative default
# so dRep does not launch checkM on 1,247 genomes.
FALLBACK_COMP = 95.0
FALLBACK_CONT = 2.0


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


refs = read_tsv(REFS)
print("=" * 74)
print("REFERENCE SET")
print("=" * 74)
print("  rows in %s: %d" % (os.path.basename(REFS), len(refs)))
print("  in host-filtered set: %s"
      % dict(Counter(r["in_host_filtered_set"] for r in refs)))
print("  Using ALL 1,247, not the 1,007 host-filtered subset. The filter is")
print("  what made the family trees STALE.")

# ------------------------------------------------- completeness metadata
qual = {}
if os.path.exists(META):
    print()
    print("  reading quality from %s" % META)
    rows = read_tsv(META)
    hdr = list(rows[0].keys()) if rows else []
    acc_c = next((c for c in ("accession", "ncbi_genbank_assembly_accession")
                  if c in hdr), None)
    comp_c = next((c for c in ("checkm_completeness", "checkm2_completeness")
                   if c in hdr), None)
    cont_c = next((c for c in ("checkm_contamination", "checkm2_contamination")
                   if c in hdr), None)
    print("  columns used: %s / %s / %s" % (acc_c, comp_c, cont_c))
    if acc_c and comp_c and cont_c:
        for r in rows:
            a = r[acc_c].strip().replace("GB_", "").replace("RS_", "")
            try:
                qual[a] = (float(r[comp_c]), float(r[cont_c]))
            except ValueError:
                pass
    print("  quality rows parsed: %d" % len(qual))
else:
    print()
    print("  NO release metadata at %s" % META)
    print("  Falling back to comp=%.1f con=%.1f for references."
          % (FALLBACK_COMP, FALLBACK_CONT))
    print("  This is a stated assumption, not measured, and must go in Methods.")

# --------------------------------------------------------------- stage
print()
print("=" * 74)
print("STAGING")
print("=" * 74)
before = len([e for e in os.scandir(GDIR) if e.name.endswith(".fa")])
print("  already staged (548 expected): %d" % before)

n_new = n_skip = n_fail = 0
staged_refs = {}
for r in refs:
    a = r["accession"].strip()
    src = r["path"].strip()
    safe = "ref__" + a
    dst = os.path.join(GDIR, safe + ".fa")
    if os.path.exists(dst) and os.path.getsize(dst) > 0:
        staged_refs[safe] = (a, src)
        n_skip += 1
        continue
    try:
        with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        if os.path.getsize(dst) == 0:
            raise IOError("wrote zero bytes")
        staged_refs[safe] = (a, src)
        n_new += 1
    except Exception as ex:
        n_fail += 1
        if n_fail <= 5:
            print("    FAIL %s: %s" % (a, ex))
        if os.path.exists(dst):
            os.remove(dst)

print("  decompressed %d, already present %d, failed %d" % (n_new, n_skip, n_fail))

after = [e for e in os.scandir(GDIR) if e.name.endswith(".fa")]
tot = sum(e.stat().st_size for e in after)
zero = [e.name for e in after if e.stat().st_size == 0]
print("  total staged files: %d, %.1f GB" % (len(after), tot / 1e9))
print("  zero-byte: %d %s" % (len(zero), zero[:5]))

# ---------------------------------------------------- append metadata
print()
print("=" * 74)
print("APPENDING TO genome_arms.tsv AND genome_info.csv")
print("=" * 74)
existing_map = set()
with open(MAP) as fh:
    fh.readline()
    for line in fh:
        existing_map.add(line.split("\t")[0])
n_map = 0
with open(MAP, "a") as f:
    for safe, (a, src) in sorted(staged_refs.items()):
        if safe + ".fa" in existing_map:
            continue
        f.write("%s\t%s\t%s\t%s\n" % (safe + ".fa", "gtdb_ref", a, src))
        n_map += 1
print("  appended %d rows to genome_arms.tsv" % n_map)

existing_info = set()
with open(INFO) as fh:
    fh.readline()
    for line in fh:
        existing_info.add(line.split(",")[0])
n_info = 0
n_default = 0
with open(INFO, "a") as f:
    for safe, (a, src) in sorted(staged_refs.items()):
        if safe + ".fa" in existing_info:
            continue
        c, x = qual.get(a, (FALLBACK_COMP, FALLBACK_CONT))
        if a not in qual:
            n_default += 1
        f.write("%s.fa,%.2f,%.2f\n" % (safe, c, x))
        n_info += 1
print("  appended %d rows to genome_info.csv (%d used the fallback)"
      % (n_info, n_default))

n_files = len([e for e in os.scandir(GDIR) if e.name.endswith(".fa")])
n_rows = sum(1 for _ in open(INFO)) - 1
print()
print("  staged .fa files : %d" % n_files)
print("  genome_info rows : %d" % n_rows)
if n_files != n_rows:
    print("  MISMATCH. dRep would run checkM on the unlisted genomes. Fix first.")
else:
    print("  Counts agree, dRep will skip checkM.")

by_arm = Counter()
with open(MAP) as fh:
    fh.readline()
    for line in fh:
        by_arm[line.split("\t")[1]] += 1
print()
print("  final composition: %s" % dict(by_arm))

# SENTINEL_END
