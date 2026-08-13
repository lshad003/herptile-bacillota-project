#!/usr/bin/env python3
# Two additional neighborhood genomes are staged for annotation.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/stage_angelakisella_emapper_gap2.py
# Output: work/angelakisella_emapper_gap2/faa/, work/angelakisella_emapper_gap2/tasks.tsv
"""
Stages the 2 additional genomes needed for the 61-tip Angelakisella
neighborhood figure: the Youngblut Angelakisella member and one
unannotated Heteroruminococcus reference.

NOTE: the Youngblut genome ID contains a double underscore
(SAMEA104404100__metabat2_low_PE.067). IDs are handled as literal
strings; nothing splits on "__" here.

Writes work/angelakisella_emapper_gap2/faa/ and tasks.tsv.
Refuses to overwrite an existing task list.
"""

import os
import shutil
import sys

BASE = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
PROTEOMES = os.path.join(BASE, "work/all_arms_pangenome/proteomes")
STAGE = os.path.join(BASE, "work/angelakisella_emapper_gap2")
FAA = os.path.join(STAGE, "faa")
TASKS = os.path.join(STAGE, "tasks.tsv")

MISSING = [
    ("youngblut", "SAMEA104404100__metabat2_low_PE.067", "eggnog"),
    ("gtdb_ref", "GCA_904387055.1", "eggnog_refs"),
]

if os.path.exists(TASKS):
    sys.exit("STOP: task list exists, refusing to overwrite:\n" + TASKS)

os.makedirs(FAA, exist_ok=True)

rows = []
for arm, gid, outdir in MISSING:
    src = os.path.join(PROTEOMES, "%s__%s.faa" % (arm, gid))
    if not os.path.isfile(src) or os.path.getsize(src) == 0:
        sys.exit("STOP: missing or empty proteome:\n" + src)
    dst = os.path.join(FAA, gid + ".faa")
    shutil.copyfile(src, dst)
    n = 0
    with open(dst) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    rows.append((gid, outdir, n))
    print("staged %-40s %-12s %6d proteins" % (gid, outdir, n))

with open(TASKS, "w") as out:
    out.write("index\tgenome\toutdir\n")
    for i, (gid, outdir, n) in enumerate(rows, 1):
        out.write("%d\t%s\t%s\n" % (i, gid, outdir))

print("")
print("task list written:", TASKS)
print("tasks: %d, array should be 1-%d" % (len(rows), len(rows)))

# STAGE_ANGELAKISELLA_EMAPPER_GAP2_V1
