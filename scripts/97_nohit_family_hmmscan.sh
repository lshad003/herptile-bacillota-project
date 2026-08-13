#!/usr/bin/env bash
# Domain content of the family determined by hmmscan against Pfam (FkbH_N).
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/angelakisella_nohit_hmmscan.sh
# Output: work/angelakisella_diamond/nohit_probe_pfam_domtbl.txt
# hmmscan of the 597 aa probe (herptile__UHM967.23060_R.bin.14|k141_274340_2)
# against Pfam 38.2 to determine the ACTUAL domain content of the
# 532-member family. The emapper Description (HAD phosphatase) and PFAMs
# (8 PKS/NRPS domains) cannot both fit a 593 aa protein; one is OG-level
# transfer. This settles it with a direct domain scan.
# Single query: login node is fine.

WORK=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/work/angelakisella_diamond
Q=$WORK/nohit_probe.faa
PFAM=/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm
OUT=$WORK/nohit_probe_pfam.txt
TBL=$WORK/nohit_probe_pfam_domtbl.txt

module load hmmer/3.3.2

HMMSCAN=$(command -v hmmscan)
if [ -z "$HMMSCAN" ]; then
    echo "STOP: hmmscan not found after module load hmmer/3.3.2" >&2
    exit 1
fi
echo "hmmscan: $HMMSCAN"

if [ ! -s "$Q" ]; then
    echo "STOP: probe FASTA missing: $Q" >&2
    exit 1
fi

if [ -s "$TBL" ]; then
    echo "STOP: output exists, refusing to overwrite: $TBL" >&2
    exit 1
fi

$HMMSCAN --cut_ga --domtblout "$TBL" -o "$OUT" "$PFAM" "$Q" || exit 1

echo "--- domains at Pfam gathering thresholds ---"
grep -v '^#' "$TBL" | awk '{printf "%-22s %-12s ali %4d-%4d  ievalue %s\n", $1, $2, $18, $19, $13}'
echo "---"
echo "Columns: Pfam domain, accession, alignment coords on the 597 aa"
echo "protein, independent e-value. --cut_ga means only domains passing"
echo "the curated Pfam thresholds appear; this is the protein's real"
echo "domain architecture, not OG-transferred annotation."
# ANGELAKISELLA_NOHIT_HMMSCAN_V1
