#!/usr/bin/env python3
# Individual Pfam models are extracted by streaming, since the installed database carries no index.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/extract_pfam_models.py
# Output: work/focal_genus_pangenome/pyrg_control/hmm/target_models.hmm
# Pulls six Pfam models out of the flat Pfam-A.hmm by streaming.
# hmmfetch cannot be used on /srv/projects/db/pfam/2021-11-25-Pfam35.0:
# it requires Pfam-A.hmm.ssi, the directory holds only the hmmpress
# binaries (.h3f .h3i .h3m .h3p), and it is owned by pkgadmin so the
# index cannot be built there. Jobs 27243541 and 27244739 both fell
# back to scanning all of Pfam 35.0 and timed out.

import os, sys

PFAM = "/srv/projects/db/pfam/2021-11-25-Pfam35.0/Pfam-A.hmm"
OUTDIR = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent/work/focal_genus_pangenome/pyrg_control/hmm"
OUT = os.path.join(OUTDIR, "target_models.hmm")

WANT = [
    ("PF00334", "NDK, nucleoside diphosphate kinase, universal control"),
    ("PF01773", "Nucleos_tra2_N, NupC nucleoside transporter"),
    ("PF07670", "Gate, NupC gate domain"),
    ("PF03825", "Nuc_H_symport, NupG nucleoside:H+ symporter"),
    ("PF06418", "CTP_synth_N, pyrG"),
    ("PF00117", "GATase, glutamine amidotransferase, positive control"),
]
WANTED = dict(WANT)

if not os.path.exists(PFAM):
    sys.exit("Pfam-A.hmm not found at %s" % PFAM)

os.makedirs(OUTDIR, exist_ok=True)

found = {}
buf = []
acc = None
name = None
n_records = 0

with open(PFAM, "r", errors="replace") as fh, open(OUT, "w") as out:
    for line in fh:
        if line.startswith("HMMER3/"):
            buf = [line]
            acc = None
            name = None
            continue
        if not buf:
            continue
        buf.append(line)
        if line.startswith("ACC "):
            acc = line.split()[1].split(".")[0]
        elif line.startswith("NAME "):
            name = line.split()[1]
        elif line.rstrip() == "//":
            n_records += 1
            if acc in WANTED and acc not in found:
                out.writelines(buf)
                found[acc] = name
            buf = []

print("records scanned : %d" % n_records)
print("models written  : %d of %d" % (len(found), len(WANT)))
for a, d in WANT:
    if a in found:
        print("  FOUND    %s  %-16s  %s" % (a, found[a], d))
    else:
        print("  MISSING  %s  %s" % (a, d))
print("output: %s" % OUT)
if len(found) != len(WANT):
    sys.exit("NOT ALL MODELS FOUND, do not submit job 47")
print("EXTRACT_PFAM_MODELS_V2_20260805_COMPLETE")
