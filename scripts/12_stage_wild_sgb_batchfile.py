# FASTA paths are resolved for the wild SGB representatives and written as a GTDB-Tk batch file.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/stage_wild_sgb_batchfile.py
# Output: data/tasks/wild_sgb_gtdbtk_batchfile.tsv, data/tasks/wild_sgb_missing_fasta.tsv
# STAGE_WILD_SGB_V1_20260803
import os
from collections import Counter

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
MANIFEST = os.path.join(ROOT, "data/sgb_manifest.tsv")
OUT_BATCH = os.path.join(ROOT, "data/tasks/wild_sgb_gtdbtk_batchfile.tsv")
OUT_MISS = os.path.join(ROOT, "data/tasks/wild_sgb_missing_fasta.tsv")

# Directories that plausibly hold herptile bin FASTAs, most likely first.
SEARCH = [
    "work/gunc_input",
    "work/drep_input",
    "results/drep_herptile_95ani_2229",
    "work/sgb_representatives",
    "data/mags",
    "work",
]
EXT = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")


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


def strip_ext(name):
    for e in (".fa.gz", ".fna.gz", ".fasta.gz", ".fa", ".fna", ".fasta"):
        if name.endswith(e):
            return name[:-len(e)]
    return name


rows = read_tsv(MANIFEST)
wild = [r for r in rows if r["has_wild"].strip().lower() == "yes"]
reps = {r["representative"].strip(): r for r in wild}
print("=" * 74)
print("TARGET SET")
print("=" * 74)
print("  wild SGBs in manifest: %d" % len(wild))
print("  distinct representatives: %d" % len(reps))
print("  by family:")
for fam, n in Counter(r["family"].strip() for r in wild).most_common(6):
    print("     %-26s %4d" % (fam, n))

# ------------------------------------------------------------ locate fastas
print()
print("=" * 74)
print("LOCATING FASTA FILES")
print("=" * 74)
found = {}
scanned = 0
for rel in SEARCH:
    base = os.path.join(ROOT, rel)
    if not os.path.isdir(base):
        print("  skip (absent): %s" % rel)
        continue
    hits_here = 0
    for dirpath, dirnames, filenames in os.walk(base):
        dirnames[:] = [d for d in dirnames if d not in ("venv", ".git", "node_modules")]
        for fn in filenames:
            if not fn.endswith(EXT):
                continue
            scanned += 1
            stem = strip_ext(fn)
            if stem in reps and stem not in found:
                found[stem] = os.path.join(dirpath, fn)
                hits_here += 1
    print("  %-44s %5d representatives matched" % (rel, hits_here))
    if len(found) == len(reps):
        print("  all representatives located, stopping search")
        break

print()
print("  fasta files scanned: %d" % scanned)
print("  representatives located: %d of %d" % (len(found), len(reps)))

missing = sorted(set(reps) - set(found))
if missing:
    print("  MISSING: %d" % len(missing))
    for m in missing[:15]:
        print("     %s" % m)
    if len(missing) > 15:
        print("     ... and %d more" % (len(missing) - 15))

# --------------------------------------------------------- size sanity check
if found:
    sizes = []
    zero = []
    for g, p in found.items():
        try:
            s = os.path.getsize(p)
        except OSError:
            s = -1
        sizes.append(s)
        if s <= 0:
            zero.append(g)
    sizes.sort()
    print()
    print("  fasta size MB: min %.2f  median %.2f  max %.2f"
          % (sizes[0] / 1e6, sizes[len(sizes) // 2] / 1e6, sizes[-1] / 1e6))
    print("  zero-byte or unreadable: %d %s" % (len(zero), zero[:5] if zero else ""))
    print("  Two zero-byte assemblies caused the 36 missing proteomes in July.")

    gz = sum(1 for p in found.values() if p.endswith(".gz"))
    print("  gzipped: %d of %d. GTDB-Tk accepts gzipped input." % (gz, len(found)))

    dirs = Counter(os.path.dirname(p) for p in found.values())
    print()
    print("  representatives came from %d directories:" % len(dirs))
    for d, n in dirs.most_common(5):
        print("     %5d  %s" % (n, d))

os.makedirs(os.path.dirname(OUT_BATCH), exist_ok=True)
with open(OUT_BATCH, "w") as f:
    for g in sorted(found):
        f.write("%s\t%s\n" % (found[g], g))
print()
print("wrote %s  (%d rows, no header, path<TAB>genome_id)" % (OUT_BATCH, len(found)))

if missing:
    with open(OUT_MISS, "w") as f:
        f.write("representative\tsgb\tfamily\n")
        for g in missing:
            r = reps[g]
            f.write("%s\t%s\t%s\n" % (g, r["sgb"], r["family"]))
    print("wrote %s  (%d rows)" % (OUT_MISS, len(missing)))

# ------------------------------------------------- youngblut herptile check
print()
print("=" * 74)
print("YOUNGBLUT HERPTILE HOSTS: ARE THEY RECOVERABLE?")
print("=" * 74)
print("  The 393 downloaded excluded herptile and fish hosts. Checking whether")
print("  the fetch table still lists them, which decides if they can ride along.")
qc = os.path.join(ROOT, "data/youngblut/youngblut_fetch_qc.tsv")
if os.path.exists(qc):
    yb = read_tsv(qc)
    print("  rows in youngblut_fetch_qc.tsv: %d" % len(yb))
    print("  distinct hosts: %d" % len({r["host"].strip() for r in yb if r["host"].strip()}))
    print("  This file holds only what was DOWNLOADED, so herptile hosts will be")
    print("  absent. Recovering them needs the original Youngblut supplementary")
    print("  table, not this file. Reporting host list for the record:")
    for h, n in Counter(r["host"].strip() for r in yb).most_common(50):
        print("     %4d  %s" % (n, h))
else:
    print("  MISSING %s" % qc)

# SENTINEL_END
