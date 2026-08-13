#!/usr/bin/env python3
# Sixteen genomes lacking eggNOG annotation are staged for annotation.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/stage_angelakisella_emapper_gap.py
# Output: work/angelakisella_emapper_gap/faa/, work/angelakisella_emapper_gap/tasks.tsv
"""
Stages the 16 Angelakisella genomes lacking eggNOG annotation (10 EHI
newt + 6 GTDB references, from angelakisella_annot_coverage.py) for an
emapper array job.

Copies each proteome from work/all_arms_pangenome/proteomes/ into
work/angelakisella_emapper_gap/faa/ under its bare genome ID, and writes
a task list mapping array index -> genome ID -> destination directory.

Refuses to overwrite an existing task list.
"""

import os
import shutil
import sys

BASE = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
PROTEOMES = os.path.join(BASE, "work/all_arms_pangenome/proteomes")
STAGE = os.path.join(BASE, "work/angelakisella_emapper_gap")
FAA = os.path.join(STAGE, "faa")
TASKS = os.path.join(STAGE, "tasks.tsv")

MISSING = [
    ("ehi_amphibian", "EHM034541", "eggnog"),
    ("ehi_amphibian", "EHM034833", "eggnog"),
    ("ehi_amphibian", "EHM035004", "eggnog"),
    ("ehi_amphibian", "EHM047839", "eggnog"),
    ("ehi_amphibian", "EHM047946", "eggnog"),
    ("ehi_amphibian", "EHM058353", "eggnog"),
    ("ehi_amphibian", "EHM058845", "eggnog"),
    ("ehi_amphibian", "EHM059119", "eggnog"),
    ("ehi_amphibian", "EHM059992", "eggnog"),
    ("ehi_amphibian", "EHM062341", "eggnog"),
    ("gtdb_ref", "GCA_900552845.1", "eggnog_refs"),
    ("gtdb_ref", "GCA_904420255.1", "eggnog_refs"),
    ("gtdb_ref", "GCA_937914895.1", "eggnog_refs"),
    ("gtdb_ref", "GCA_949285935.1", "eggnog_refs"),
    ("gtdb_ref", "GCA_949298535.1", "eggnog_refs"),
    ("gtdb_ref", "GCF_900104675.1", "eggnog_refs"),
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
    print("staged %-18s %-12s %6d proteins" % (gid, outdir, n))

with open(TASKS, "w") as out:
    out.write("index\tgenome\toutdir\n")
    for i, (gid, outdir, n) in enumerate(rows, 1):
        out.write("%d\t%s\t%s\n" % (i, gid, outdir))

print("")
print("task list written:", TASKS)
print("tasks: %d, array should be 1-%d" % (len(rows), len(rows)))

# STAGE_ANGELAKISELLA_EMAPPER_GAP_V1
