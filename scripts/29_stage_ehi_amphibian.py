# EHI newt genomes are staged as the independent amphibian comparison arm.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/stage_ehi_amphibian.py
# Output: results/ehi_amphibian_manifest.tsv, data/tasks/ehi_amphibian_gtdbtk_batchfile.tsv
# STAGE_EHI_AMPHIBIAN_V1_20260804
import os
from collections import Counter

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
LIST = os.path.join(CH3, "data/ehi_amphibian_genome_list.tsv")
FADIR = os.path.join(CH3, "data/ehi_2025/mags/amphibian_fa")

WORK = os.path.join(ROOT, "work/ehi_amphibian")
LINKS = os.path.join(WORK, "genomes")
OUT_MAN = os.path.join(ROOT, "results/ehi_amphibian_manifest.tsv")
OUT_BATCH = os.path.join(ROOT, "data/tasks/ehi_amphibian_gtdbtk_batchfile.tsv")

PHYLUM = "p__Bacillota_A"
COMP_MIN = 70.0
CONT_MAX = 10.0


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


def tofloat(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


rows = read_tsv(LIST)
print("=" * 74)
print("EHI AMPHIBIAN GENOME LIST")
print("=" * 74)
print("  rows: %d" % len(rows))
print("  hosts: %s" % dict(Counter(r["host_species"].strip() for r in rows)))
print()
print("  by phylum (top 8):")
for p, n in Counter(r["phylum"].strip() for r in rows).most_common(8):
    print("    %5d  %s" % (n, p))

target = [r for r in rows if r["phylum"].strip() == PHYLUM]
print()
print("  %s: %d" % (PHYLUM, len(target)))
print("  by host: %s" % dict(Counter(r["host_species"].strip() for r in target)))

# quality filter matched to the pooled dRep settings, so the arm is comparable
kept = []
dropped = 0
for r in target:
    c, x = tofloat(r["completeness"]), tofloat(r["contamination"])
    if c is None or x is None:
        dropped += 1
        continue
    if c >= COMP_MIN and x <= CONT_MAX:
        kept.append(r)
    else:
        dropped += 1
print()
print("  quality filter comp>=%.0f con<=%.0f (matching the pooled dRep run)"
      % (COMP_MIN, CONT_MAX))
print("    kept %d, dropped %d" % (len(kept), dropped))
if kept:
    cs = sorted(tofloat(r["completeness"]) for r in kept)
    print("    completeness: min %.1f median %.1f max %.1f"
          % (cs[0], cs[len(cs) // 2], cs[-1]))

# ------------------------------------------------------------- locate fastas
print()
print("=" * 74)
print("LOCATING FASTA FILES")
print("=" * 74)
if not os.path.isdir(FADIR):
    raise SystemExit("missing %s" % FADIR)
present = {}
for e in os.scandir(FADIR):
    if e.is_file() and e.name.endswith((".fa.gz", ".fa", ".fna.gz", ".fna")):
        stem = e.name.split(".")[0]
        present[stem] = e.path
print("  files in amphibian_fa: %d" % len(present))

found, missing = [], []
for r in kept:
    g = r["genome_id"].strip()
    if g in present:
        found.append((r, present[g]))
    else:
        missing.append(g)
print("  located: %d, missing: %d %s" % (len(found), len(missing), missing[:5]))

zero = [p for _, p in found if os.path.getsize(p) == 0]
print("  zero-byte: %d" % len(zero))

for d in (WORK, LINKS, os.path.dirname(OUT_BATCH)):
    if not os.path.isdir(d):
        os.makedirs(d)

n_link = 0
for r, p in found:
    g = r["genome_id"].strip()
    dst = os.path.join(LINKS, g + ".fa.gz")
    if not os.path.exists(dst):
        os.symlink(p, dst)
        n_link += 1
print("  symlinks created: %d in %s" % (n_link, LINKS))

with open(OUT_MAN, "w") as f:
    f.write("genome_id\tassembly_id\thost_species\tphylum\t"
            "completeness\tcontamination\tfasta\n")
    for r, p in found:
        f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                % (r["genome_id"].strip(), r["assembly_id"].strip(),
                   r["host_species"].strip(), r["phylum"].strip(),
                   r["completeness"].strip(), r["contamination"].strip(), p))
print("  wrote %s (%d rows)" % (OUT_MAN, len(found)))

with open(OUT_BATCH, "w") as f:
    for r, p in found:
        f.write("%s\t%s\n" % (p, r["genome_id"].strip()))
print("  wrote %s (%d rows)" % (OUT_BATCH, len(found)))

print()
print("=" * 74)
print("WHY THIS ARM MATTERS")
print("=" * 74)
print("  Step 3 currently makes an ABSENCE claim: 15 genera holding 97 herptile")
print("  SGBs are in GTDB but not recovered from any endotherm gut.")
print("  EHI amphibians are newts (Caudata), the same order as the salamanders")
print("  that make up 52%% of the wild catalog, from a different lab and")
print("  different sites. If those genera appear there and stay absent from")
print("  EHI's own mammals, the claim becomes a PRESENCE claim, which needs")
print("  far less hedging.")
print()
print("  Next: GTDB-Tk classify these to find the Ruminococcaceae, then add")
print("  them to the pooled dRep as a fifth arm.")

# SENTINEL_END
