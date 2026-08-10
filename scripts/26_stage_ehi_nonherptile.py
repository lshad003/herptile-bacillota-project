#!/usr/bin/env python3
# Non-herptile EHI Bacillota_A genomes are staged as the mammal comparison arm.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/stage_ehi_nonherptile.py
# Output: work/ehi_nonherptile/genomes/, results/ehi_nonherptile_manifest.tsv
import os, sys, gzip
from collections import Counter

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
TAB = B + "/ehi_validation/results/ehi_bacillota_a_hq_genomes.tsv"
MAGD = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution/data/ehi_2025/mags"
SUBS = ("nonherptile_fa", "nonamph_fa", "amphibian_fa")
STAGE = B + "/work/ehi_nonherptile/genomes"
OUT = B + "/results/ehi_nonherptile_manifest.tsv"
BATCH = B + "/data/tasks/ehi_gtdbtk_batchfile.tsv"

KEEP = ("Mammalia", "Aves")

if not os.path.exists(TAB):
    print("MISSING:", TAB); sys.exit(1)

rows = []
with open(TAB) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    I = {k: h.index(k) for k in
         ("genome_id", "family", "host_species", "host_class", "host_order",
          "completeness", "contamination", "size_mb", "contigs")}
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(I.values()) or p[I["host_class"]] not in KEEP:
            continue
        rows.append(dict(g=p[I["genome_id"]], fam=p[I["family"]],
                         sp=p[I["host_species"]], cl=p[I["host_class"]],
                         od=p[I["host_order"]],
                         comp=float(p[I["completeness"]]),
                         cont=float(p[I["contamination"]]),
                         size=float(p[I["size_mb"]]),
                         ctg=int(p[I["contigs"]])))
print("non-herptile Bacillota_A: %d" % len(rows))
print("by class:", dict(Counter(r["cl"] for r in rows)))

index = {}
for s in SUBS:
    d = os.path.join(MAGD, s)
    if not os.path.isdir(d):
        print("  absent %s" % d); continue
    n = 0
    for e in os.scandir(d):
        if e.name.startswith("._"):
            continue
        if e.name.endswith((".fa.gz", ".fa", ".fna.gz", ".fna")):
            k = e.name.split(".")[0]
            if k not in index:
                index[k] = e.path; n += 1
    print("  %-16s indexed %d" % (s, n))
print("  total indexed: %d" % len(index))

have = [r for r in rows if r["g"] in index]
miss = [r for r in rows if r["g"] not in index]
print()
print("  found %d | missing %d" % (len(have), len(miss)))
for r in miss[:5]:
    print("     MISS %s" % r["g"])
if not have:
    print("STOP: nothing staged."); sys.exit(1)

if not os.path.isdir(STAGE):
    os.makedirs(STAGE)

ok = []; bad = []
for r in have:
    dst = os.path.join(STAGE, r["g"] + ".fa.gz")
    if not os.path.exists(dst):
        try:
            os.symlink(index[r["g"]], dst)
        except OSError as e:
            bad.append((r["g"], type(e).__name__)); continue
    try:
        with gzip.open(dst, "rt") as fh:
            first = fh.readline()
            body = len(fh.read(4096))
        if not first.startswith(">") or body < 500:
            bad.append((r["g"], "bad FASTA")); continue
    except Exception as e:
        bad.append((r["g"], "gzip " + type(e).__name__)); continue
    ok.append(r)
print("  staged and readable: %d" % len(ok))
if bad:
    print("  UNUSABLE %d:" % len(bad))
    for g, w in bad[:8]:
        print("     %s %s" % (g, w))

with open(OUT, "w") as f:
    f.write("genome_id\tehi_family\thost_species\thost_class\thost_order\t"
            "completeness\tcontamination\tsize_mb\tcontigs\tfasta\n")
    for r in ok:
        f.write("%s\t%s\t%s\t%s\t%s\t%.2f\t%.3f\t%.4f\t%d\t%s\n"
                % (r["g"], r["fam"], r["sp"], r["cl"], r["od"], r["comp"],
                   r["cont"], r["size"], r["ctg"],
                   os.path.join(STAGE, r["g"] + ".fa.gz")))
td = os.path.dirname(BATCH)
if not os.path.isdir(td):
    os.makedirs(td)
with open(BATCH, "w") as f:
    for r in ok:
        f.write("%s\t%s\n" % (os.path.join(STAGE, r["g"] + ".fa.gz"), r["g"]))

print()
print("wrote", OUT)
print("wrote", BATCH, "(GTDB-Tk batchfile, %d genomes)" % len(ok))
print()
print("EHI family labels (GTDB-Tk 2.3.0 / r214), top 8:")
for fam, c in Counter(r["fam"] for r in ok).most_common(8):
    print("   %-26s %d" % (fam, c))
print()
print("Ruminococcaceae staged: %d"
      % sum(1 for r in ok if r["fam"] == "f__Ruminococcaceae"))
print("completeness mean %.1f | contamination mean %.2f"
      % (sum(r["comp"] for r in ok) / len(ok),
         sum(r["cont"] for r in ok) / len(ok)))
print()
print("These labels are EHI's r214 call. GTDB-Tk r226 will reassign them,")
print("and the r226 genus is what the analysis uses.")
print("DONE_STAGE_EHI")
# SENTINEL_END
