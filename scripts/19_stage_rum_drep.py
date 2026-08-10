# Genomes are staged into the five-arm Ruminococcaceae dereplication input set.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/stage_rum_drep.py
# Output: work/pooled_drep_rum/genomes/, genome_arms.tsv, genome_info.csv
# STAGE_RUM_DREP_V1_20260803
import os, gzip, shutil
from collections import Counter

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
WORK = os.path.join(ROOT, "work/pooled_drep_rum")
GDIR = os.path.join(WORK, "genomes")
OUT_MAP = os.path.join(WORK, "genome_arms.tsv")
OUT_INFO = os.path.join(WORK, "genome_info.csv")

SGB_MANIFEST = os.path.join(ROOT, "data/sgb_manifest.tsv")
EHI_MANIFEST = os.path.join(ROOT, "results/ehi_nonherptile_manifest.tsv")
EHI_SUM = os.path.join(ROOT, "results/gtdbtk_ehi_r220_classify/gtdbtk.bac120.summary.tsv")
YB_SUM = os.path.join(ROOT, "results/gtdbtk_youngblut_r220/gtdbtk.bac120.summary.tsv")
YB_BATCH = os.path.join(ROOT, "data/youngblut/youngblut_gtdbtk_batchfile.tsv")
YB_QC = os.path.join(ROOT, "data/youngblut/youngblut_fetch_qc.tsv")
HERP_DIR = os.path.join(ROOT, "results/drep_herptile_95ani_2229/dereplicated_genomes")

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


for d in (WORK, GDIR):
    if not os.path.isdir(d):
        os.makedirs(d)

targets = {}   # safe_name -> (source_path, arm, genome_id, comp, cont)

# ------------------------------------------------------------- herptile
wild = [r for r in read_tsv(SGB_MANIFEST)
        if r["has_wild"].strip().lower() == "yes"
        and r["family"].strip() == FAM]
print("=" * 74)
print("BUILDING POOLED %s SET" % FAM.upper())
print("=" * 74)
print("  herptile wild %s SGBs: %d" % (FAM, len(wild)))
for r in wild:
    rep = r["representative"].strip()
    for e in (".fa", ".fna", ".fasta"):
        p = os.path.join(HERP_DIR, rep + e)
        if os.path.exists(p):
            targets["herp__" + rep] = (p, "herptile", rep,
                                       r["rep_completeness"], r["rep_contamination"])
            break

# ------------------------------------------------------------- ehi
ehi_fam = {}
for r in read_tsv(EHI_SUM):
    ehi_fam[r["user_genome"].strip()] = parse_tax(r["classification"]).get("f", "")
ehi_rows = [r for r in read_tsv(EHI_MANIFEST)
            if ehi_fam.get(r["genome_id"].strip(), "") == FAM]
print("  ehi %s MAGs: %d" % (FAM, len(ehi_rows)))
for r in ehi_rows:
    g = r["genome_id"].strip()
    p = r.get("fasta", "").strip()
    if p and os.path.exists(p):
        targets["ehi__" + g] = (p, "ehi", g, r["completeness"], r["contamination"])

# ------------------------------------------------------------- youngblut
yb_fam = {}
for r in read_tsv(YB_SUM):
    yb_fam[r["user_genome"].strip()] = parse_tax(r["classification"]).get("f", "")
yb_qc = {r["genome"].strip(): r for r in read_tsv(YB_QC)}
yb_path = {}
with open(YB_BATCH) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            yb_path[f[1].strip()] = f[0].strip()
yb_ids = [g for g, fam in yb_fam.items() if fam == FAM]
print("  youngblut %s SGBs: %d" % (FAM, len(yb_ids)))
for g in yb_ids:
    p = yb_path.get(g, "")
    if p and os.path.exists(p):
        q = yb_qc.get(g, {})
        targets["yb__" + g] = (p, "youngblut", g,
                               q.get("completeness", ""), q.get("contamination", ""))

print()
c = Counter(v[1] for v in targets.values())
print("  located: %s  total %d" % (dict(c), len(targets)))
missing_h = len(wild) - c.get("herptile", 0)
missing_e = len(ehi_rows) - c.get("ehi", 0)
missing_y = len(yb_ids) - c.get("youngblut", 0)
if missing_h or missing_e or missing_y:
    print("  MISSING herptile %d, ehi %d, youngblut %d" % (missing_h, missing_e, missing_y))

# --------------------------------------------------- stage uncompressed
print()
print("=" * 74)
print("STAGING UNCOMPRESSED COPIES (dRep cannot read .gz)")
print("=" * 74)
print("  destination: %s" % GDIR)
n_gz = n_cp = n_skip = 0
for safe, (src, arm, gid, comp, cont) in sorted(targets.items()):
    dst = os.path.join(GDIR, safe + ".fa")
    if os.path.exists(dst) and os.path.getsize(dst) > 0:
        n_skip += 1
        continue
    if src.endswith(".gz"):
        with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        n_gz += 1
    else:
        shutil.copyfile(src, dst)
        n_cp += 1
print("  decompressed %d, copied %d, already present %d" % (n_gz, n_cp, n_skip))

staged = [e for e in os.scandir(GDIR) if e.is_file() and e.name.endswith(".fa")]
tot = sum(e.stat().st_size for e in staged)
zero = [e.name for e in staged if e.stat().st_size == 0]
print("  staged files: %d, total %.1f MB" % (len(staged), tot / 1e6))
print("  zero-byte: %d %s" % (len(zero), zero[:5]))

# ------------------------------------------------------------- metadata
with open(OUT_MAP, "w") as f:
    f.write("staged_name\tarm\tgenome_id\tsource_path\n")
    for safe, (src, arm, gid, comp, cont) in sorted(targets.items()):
        f.write("%s\t%s\t%s\t%s\n" % (safe + ".fa", arm, gid, src))

# dRep genomeInfo: uses OUR completeness/contamination so it does not rerun checkM
n_info = 0
with open(OUT_INFO, "w") as f:
    f.write("genome,completeness,contamination\n")
    for safe, (src, arm, gid, comp, cont) in sorted(targets.items()):
        try:
            cv = float(comp); xv = float(cont)
        except ValueError:
            continue
        f.write("%s.fa,%.2f,%.2f\n" % (safe, cv, xv))
        n_info += 1
print()
print("  wrote %s" % OUT_MAP)
print("  wrote %s  (%d rows)" % (OUT_INFO, n_info))
if n_info < len(targets):
    print("  WARNING: %d genomes lack completeness/contamination and would force"
          % (len(targets) - n_info))
    print("  dRep to run checkM itself. Fix before submitting.")

# SENTINEL_END
