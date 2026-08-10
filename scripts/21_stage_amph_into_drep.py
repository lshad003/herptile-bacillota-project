# Amphibian-derived genomes are added to the dereplication input set.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/stage_amph_into_drep.py
# Output: work/pooled_drep_rum/genomes/, genome_arms.tsv, genome_info.csv
# STAGE_AMPH_INTO_DREP_V1_20260804
import os, gzip, shutil
from collections import Counter

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
WORK = os.path.join(ROOT, "work/pooled_drep_rum")
GDIR = os.path.join(WORK, "genomes")
MAP = os.path.join(WORK, "genome_arms.tsv")
INFO = os.path.join(WORK, "genome_info.csv")

AMPH_SUM = os.path.join(ROOT, "results/gtdbtk_ehi_amphibian_r220/gtdbtk.bac120.summary.tsv")
AMPH_MAN = os.path.join(ROOT, "results/ehi_amphibian_manifest.tsv")

FAM = "Ruminococcaceae"


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


fam_of = {}
for r in read_tsv(AMPH_SUM):
    fam_of[r["user_genome"].strip()] = parse_tax(r["classification"]).get("f", "")
print("=" * 74)
print("EHI AMPHIBIAN CLASSIFICATION")
print("=" * 74)
print("  genomes classified: %d" % len(fam_of))
print("  top families:")
for f, n in Counter(fam_of.values()).most_common(8):
    print("    %4d  %s" % (n, f))

man = {r["genome_id"].strip(): r for r in read_tsv(AMPH_MAN)}
targets = {g: man[g] for g, f in fam_of.items() if f == FAM and g in man}
print()
print("  %s with a manifest row: %d" % (FAM, len(targets)))
print("  by host: %s" % dict(Counter(r["host_species"].strip()
                                     for r in targets.values())))
comps = sorted(float(r["completeness"]) for r in targets.values())
print("  completeness: min %.1f median %.1f max %.1f"
      % (comps[0], comps[len(comps) // 2], comps[-1]))

# --------------------------------------------------------------- stage
print()
print("=" * 74)
print("STAGING")
print("=" * 74)
before = len([e for e in os.scandir(GDIR) if e.name.endswith(".fa")])
print("  already staged: %d" % before)

n_new = n_skip = n_fail = 0
staged = {}
for g, r in sorted(targets.items()):
    src = r["fasta"].strip()
    safe = "amph__" + g
    dst = os.path.join(GDIR, safe + ".fa")
    if os.path.exists(dst) and os.path.getsize(dst) > 0:
        staged[safe] = r
        n_skip += 1
        continue
    try:
        if src.endswith(".gz"):
            with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
                shutil.copyfileobj(fi, fo)
        else:
            shutil.copyfile(src, dst)
        if os.path.getsize(dst) == 0:
            raise IOError("zero bytes written")
        staged[safe] = r
        n_new += 1
    except Exception as ex:
        n_fail += 1
        if n_fail <= 5:
            print("    FAIL %s: %s" % (g, ex))
        if os.path.exists(dst):
            os.remove(dst)
print("  staged %d, already present %d, failed %d" % (n_new, n_skip, n_fail))

# --------------------------------------------------------- append rows
existing_map = set()
with open(MAP) as fh:
    fh.readline()
    for line in fh:
        existing_map.add(line.split("\t")[0])
n_map = 0
with open(MAP, "a") as f:
    for safe, r in sorted(staged.items()):
        if safe + ".fa" in existing_map:
            continue
        f.write("%s\tehi_amphibian\t%s\t%s\n"
                % (safe + ".fa", r["genome_id"].strip(), r["fasta"].strip()))
        n_map += 1

existing_info = set()
with open(INFO) as fh:
    fh.readline()
    for line in fh:
        existing_info.add(line.split(",")[0])
n_info = 0
with open(INFO, "a") as f:
    for safe, r in sorted(staged.items()):
        if safe + ".fa" in existing_info:
            continue
        f.write("%s.fa,%.2f,%.2f\n"
                % (safe, float(r["completeness"]), float(r["contamination"])))
        n_info += 1
print("  appended %d rows to genome_arms.tsv, %d to genome_info.csv"
      % (n_map, n_info))

n_files = len([e for e in os.scandir(GDIR) if e.name.endswith(".fa")])
n_rows = sum(1 for _ in open(INFO)) - 1
print()
print("  staged .fa files : %d" % n_files)
print("  genome_info rows : %d" % n_rows)
print("  %s" % ("counts agree, dRep will skip checkM" if n_files == n_rows
                else "MISMATCH: dRep would run checkM. Fix before submitting."))

by_arm = Counter()
with open(MAP) as fh:
    fh.readline()
    for line in fh:
        by_arm[line.split("\t")[1]] += 1
print()
print("  five-arm composition: %s" % dict(by_arm))
print()
print("  The dRep output directory must be REMOVED or renamed before rerunning,")
print("  or dRep will resume from the old four-arm run.")

# SENTINEL_END
