# Proteomes are staged for carbohydrate-active enzyme annotation.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/stage_cazy_targets.py
# Output: work/cazy_focal/cazy_targets.tsv
import os, sys, csv
from collections import Counter

BASE  = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
FOCAL = os.path.join(BASE, "work/focal_genus_pangenome/proteomes")
ALLP  = os.path.join(BASE, "work/all_arms_pangenome/proteomes")
AMANI = os.path.join(BASE, "work/all_arms_pangenome/all_arms_manifest.tsv")
COH   = os.path.join(BASE, "results/unassigned_clade_coherence.tsv")
OUTD  = os.path.join(BASE, "work/cazy_focal")
LIST  = os.path.join(OUTD, "cazy_targets.tsv")

def die(m):
    print("")
    print("!" * 72)
    print("FAILED: " + m)
    print("!" * 72)
    sys.exit(1)

print("=" * 72)
print("STEP 1  THE 125 FOCAL GENOMES")
print("=" * 72)
if not os.path.isdir(FOCAL):
    die("missing " + FOCAL)
focal = {}
for e in os.scandir(FOCAL):
    if not e.name.endswith(".faa"):
        continue
    stem = e.name[:-4]
    parts = stem.split("__")
    if len(parts) < 3:
        die("unexpected focal proteome name: " + e.name)
    genus, arm, gid = parts[0], parts[1], "__".join(parts[2:])
    focal[gid] = (genus, arm)
print("  focal proteomes: %d" % len(focal))
print("  by genus: %s" % dict(Counter(v[0] for v in focal.values())))
print("  by arm  : %s" % dict(Counter(v[1] for v in focal.values())))

print("")
print("=" * 72)
print("STEP 2  THE UNASSIGNED CLADE GENOMES")
print("=" * 72)
if not os.path.exists(COH):
    die("missing " + COH)
clade = {}
for r in csv.DictReader(open(COH), delimiter="\t"):
    c = (r["clade"] or "").strip()
    for g in (r["genomes"] or "").split(";"):
        g = g.strip()
        if g:
            clade[g] = c
print("  genomes in unassigned clades: %d" % len(clade))
print("  by clade: %s" % dict(Counter(clade.values())))

print("")
print("=" * 72)
print("STEP 3  MAP ONTO THE ALL-ARMS PROTEOMES")
print("=" * 72)
if not os.path.isdir(ALLP):
    die("missing " + ALLP)
rows = list(csv.DictReader(open(AMANI), delimiter="\t"))
byid = {}
for r in rows:
    fn = "%s__%s.faa" % (r["arm"], r["genome_id"])
    p = os.path.join(ALLP, fn)
    if os.path.exists(p):
        byid[r["genome_id"]] = (p, r["arm"], r["genus"] or "unassigned")
print("  all-arms proteomes indexed: %d" % len(byid))

targets, miss_f, miss_c = [], [], []
for gid, (genus, arm) in sorted(focal.items()):
    if gid in byid:
        p, a2, g2 = byid[gid]
        targets.append((gid, p, "focal125", genus, arm, ""))
    else:
        miss_f.append(gid)
for gid, c in sorted(clade.items()):
    if gid in byid:
        p, a2, g2 = byid[gid]
        targets.append((gid, p, "unassigned_clade", g2, a2, c))
    else:
        miss_c.append(gid)

print("  focal 125 mapped     : %d of %d" % (len(focal) - len(miss_f), len(focal)))
print("  unassigned mapped    : %d of %d" % (len(clade) - len(miss_c), len(clade)))
if miss_f:
    print("  FOCAL NOT FOUND: %s" % miss_f[:6])
if miss_c:
    print("  CLADE NOT FOUND: %s" % miss_c[:6])

seen, dedup = set(), []
for t in targets:
    if t[0] in seen:
        print("  NOTE: %s is in both sets, scanned once, labelled focal125" % t[0])
        continue
    seen.add(t[0])
    dedup.append(t)
print("  unique genomes to scan: %d" % len(dedup))
print("  by set: %s" % dict(Counter(t[2] for t in dedup)))

if miss_f or miss_c:
    die("not every target maps to a proteome, refusing to scan a partial set")

if not os.path.isdir(OUTD):
    os.makedirs(OUTD)
if os.path.exists(LIST):
    print("")
    print("  NOT overwriting existing " + LIST)
else:
    with open(LIST, "w") as fh:
        fh.write("index\tgenome_id\tproteome\tset\tgenus\tarm\tclade\n")
        for i, t in enumerate(dedup):
            fh.write("%d\t%s\t%s\t%s\t%s\t%s\t%s\n" % (i, t[0], t[1], t[2], t[3], t[4], t[5]))
    print("")
    print("  wrote %s (%d rows)" % (LIST, len(dedup)))
print("")
print("  All proteomes come from the SAME Prodigal run (tonight's all-arms")
print("  array), so gene calling is uniform across both sets. The August 4")
print("  focal proteomes were NOT used.")
print("  Array size = %d tasks." % len(dedup))
print("")
print("STAGE_CAZY_TARGETS_V1_20260806 COMPLETE")
# STAGE_CAZY_TARGETS_V1_20260806
