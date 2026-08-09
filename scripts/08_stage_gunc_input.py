#!/usr/bin/env python3
# Preparation of protein calls for chimerism screening.
#
# Source: ruminococcaceae-agent/scripts/stage_gunc_input.py
# Writes: work/gunc_input/
#
# GUNC is run on protein calls rather than nucleotide sequence, so Prodigal
# is run first with -p meta. The same Prodigal version is used throughout the
# project so that gene calls are comparable between analyses.
import os, sys

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
SGB = B + "/data/sgb_manifest.tsv"
PROT = B + "/results/prodigal"
STAGE = B + "/work/gunc_input"

if not os.path.exists(SGB):
    print("MISSING:", SGB); sys.exit(1)

reps = []
with open(SGB) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    ri, fi = h.index("representative"), h.index("family")
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) > max(ri, fi):
            reps.append((p[ri], p[fi]))
print("SGB representatives: %d" % len(reps))

if not os.path.isdir(STAGE):
    os.makedirs(STAGE)

ok = 0; miss = []
for g, fam in reps:
    src = os.path.join(PROT, g + ".faa")
    if not os.path.exists(src):
        miss.append(g); continue
    dst = os.path.join(STAGE, g + ".faa")
    if not os.path.exists(dst):
        try:
            os.symlink(src, dst)
        except OSError as e:
            miss.append(g); continue
    ok += 1

print("staged: %d | missing proteome: %d" % (ok, len(miss)))
for g in miss[:10]:
    print("   MISS", g)
print("staged in", STAGE)
if miss:
    print("STOP: some representatives have no .faa. Fix before running.")
    sys.exit(1)
print("DONE_STAGE_GUNC")
# SENTINEL_END
