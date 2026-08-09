#!/usr/bin/env python3
# Recounting of catalog totals with and without the laboratory-reared arm.
#
# Source: ruminococcaceae-agent/scripts/recount_nowf.py
# Reads:  data/herptile_bacillota_A_HQ_manifest_with_source.tsv
#         data/sgb_manifest.tsv, data/sgb_manifest_nowf.tsv
#         the NCBI submission sheet, for sequencing effort
# Writes: results/recount_nowf.tsv
#
# Library and animal counts derived from the MAG manifest describe recovery
# rather than sampling effort, since the manifest contains only libraries that
# yielded at least one genome. Effort denominators are taken from the NCBI
# submission sheet, which records every library sequenced.
import os, sys
from collections import Counter, defaultdict

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
MAN = B + "/data/herptile_bacillota_A_HQ_manifest_with_source.tsv"
SGB_OLD = B + "/data/sgb_manifest.tsv"
SGB_NEW = B + "/data/sgb_manifest_nowf.tsv"
NCBI = "/bigdata/stajichlab/lshad003/ncbi-deposit/results/UHM_Metagenomes_FILLED.tsv"
OUT = B + "/results/recount_nowf.tsv"

DROP = "VIVARIUM"

def die(m):
    sys.stderr.write("FATAL: " + m + "\n"); sys.exit(1)

if os.path.exists(OUT):
    die("output exists, refusing to overwrite: " + OUT)
for p in (MAN, SGB_OLD, SGB_NEW):
    if not os.path.isfile(p):
        die("missing " + p)

def read_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return head, [dict(zip(head, l.rstrip("\n").split("\t")))
                      for l in fh if l.strip()]

mh, mrows = read_tsv(MAN)
if len(mrows) != 2229:
    die("expected 2229 manifest rows, found %d" % len(mrows))

keep = [r for r in mrows if r["source"] != DROP]
gone = [r for r in mrows if r["source"] == DROP]
if len(gone) != 64:
    die("expected 64 %s MAGs, found %d" % (DROP, len(gone)))

def summarise(rows, label):
    libs = set(r["sample_id_full"] for r in rows)
    anim = set(r["sample_id_base"] for r in rows)
    taxa = set(r["host_taxon"] for r in rows if r["host_taxon"])
    comp = [float(r["completeness"]) for r in rows]
    cont = [float(r["contamination"]) for r in rows]
    nc = sum(1 for r in rows
             if float(r["completeness"]) >= 90.0
             and float(r["contamination"]) <= 5.0)
    print("")
    print("=== %s ===" % label)
    print("  MAGs                     %d" % len(rows))
    print("  libraries yielding a MAG %d" % len(libs))
    print("  animals yielding a MAG   %d" % len(anim))
    print("  host taxa                %d" % len(taxa))
    print("  mean completeness        %.2f" % (sum(comp) / len(comp)))
    print("  mean contamination       %.2f" % (sum(cont) / len(cont)))
    print("  near-complete (>=90,<=5) %d (%.1f%%)"
          % (nc, 100.0 * nc / len(rows)))
    return dict(mags=len(rows), libs=len(libs), animals=len(anim),
                taxa=len(taxa), comp=sum(comp) / len(comp),
                cont=sum(cont) / len(cont), nc=nc)

A = summarise(mrows, "WITH THE WOOD FROG ARM")
Bk = summarise(keep, "WITHOUT THE WOOD FROG ARM")

print("")
print("=== HOST TAXA LOST ENTIRELY ===")
t_all = set(r["host_taxon"] for r in mrows if r["host_taxon"])
t_keep = set(r["host_taxon"] for r in keep if r["host_taxon"])
lost_t = sorted(t_all - t_keep)
if lost_t:
    for t in lost_t:
        n = sum(1 for r in gone if r["host_taxon"] == t)
        print("  %-40s %d MAGs" % (t, n))
else:
    print("  none")

print("")
print("=== ANIMAL TYPES, MAG COUNTS ===")
ta = Counter(r["animal_type"] for r in mrows)
tk = Counter(r["animal_type"] for r in keep)
print("  %-16s %8s %8s %8s" % ("animal_type", "with", "without", "delta"))
for k in sorted(ta, key=lambda x: -ta[x]):
    print("  %-16s %8d %8d %8d" % (k, ta[k], tk.get(k, 0), tk.get(k, 0) - ta[k]))

print("")
print("=== SOURCE SPLIT, RETAINED MAGs ===")
for k, v in sorted(Counter(r["source"] for r in keep).items(),
                   key=lambda x: -x[1]):
    print("  %-12s %d" % (k, v))

oh, orows = read_tsv(SGB_OLD)
nh, nrows = read_tsv(SGB_NEW)

def fam_table(sgbrows, magrows, famkey):
    mags = Counter(r["family"] for r in magrows)
    sgbs = Counter(r[famkey] for r in sgbrows)
    wild = Counter(r[famkey] for r in sgbrows if r["has_wild"] == "yes")
    return mags, sgbs, wild

ma, sa, wa = fam_table(orows, mrows, "family")
mb, sb, wb = fam_table(nrows, keep, "family")

print("")
print("=== BY FAMILY: MAGs / SGBs / wild SGBs ===")
print("  %-24s %18s %18s" % ("family", "with", "without"))
fams = sorted(set(sa) | set(sb), key=lambda f: -sa.get(f, 0))
for f in fams[:12]:
    print("  %-24s %6d %5d %5d %6d %5d %5d"
          % (f, ma.get(f, 0), sa.get(f, 0), wa.get(f, 0),
             mb.get(f, 0), sb.get(f, 0), wb.get(f, 0)))

print("")
print("=== SEQUENCING EFFORT, FROM THE NCBI SUBMISSION SHEET ===")
if os.path.isfile(NCBI):
    with open(NCBI) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ib = h.index("host_metagenome_bucket")
        isn = h.index("*sample_name")
        libs = defaultdict(set); anims = defaultdict(set)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(ib, isn):
                continue
            b = p[ib].strip()
            grp = "wood frog" if b.startswith("WF") else b
            libs[grp].add(p[isn].strip())
            anims[grp].add(p[isn].split(".")[0].strip())
    print("  %-14s %10s %10s" % ("bucket", "libraries", "animals"))
    for k in sorted(libs, key=lambda x: -len(libs[x])):
        print("  %-14s %10d %10d" % (k, len(libs[k]), len(anims[k])))
    tot_l = sum(len(v) for v in libs.values())
    tot_a = len(set().union(*anims.values()))
    print("  %-14s %10d %10d" % ("TOTAL", tot_l, tot_a))
    ret_l = sum(len(v) for k, v in libs.items() if k != "wood frog")
    ret_a = len(set().union(*[v for k, v in anims.items() if k != "wood frog"]))
    print("  %-14s %10d %10d" % ("RETAINED", ret_l, ret_a))
else:
    print("  NCBI sheet not readable from here; effort denominators unknown")

with open(OUT, "w") as f:
    f.write("quantity\twith_wf\twithout_wf\n")
    for k in ("mags", "libs", "animals", "taxa", "nc"):
        f.write("%s\t%s\t%s\n" % (k, A[k], Bk[k]))
    f.write("mean_completeness\t%.2f\t%.2f\n" % (A["comp"], Bk["comp"]))
    f.write("mean_contamination\t%.2f\t%.2f\n" % (A["cont"], Bk["cont"]))
    f.write("sgbs\t%d\t%d\n" % (len(orows), len(nrows)))
    f.write("wild_sgbs\t%d\t%d\n"
            % (sum(1 for r in orows if r["has_wild"] == "yes"),
               sum(1 for r in nrows if r["has_wild"] == "yes")))
print("")
print("WROTE: " + OUT)
# RECOUNT_NOWF_V1_20260808
