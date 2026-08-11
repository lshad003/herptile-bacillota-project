#!/bin/bash
# Focal proteomes are annotated with Pfam gathering cutoffs and KofamScan adaptive thresholds.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/66_annot_154_v2.sh
# Output: work/annot_154/pfam/, work/annot_154/kofam/
#SBATCH --job-name=annot154
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --array=0-153%30
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/annot154v2_%A_%a.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/annot154v2_%A_%a.err

KOFAM=/opt/linux/rocky/8.x/x86_64/pkgs/kofamscan/1.3.0/env/bin/exec_annotation
HMMSEARCH=/opt/linux/rocky/8.x/x86_64/pkgs/hmmer/3.4/bin/hmmsearch
KOPROF=/srv/projects/db/KEGG/97.0/profiles
KOLIST=/srv/projects/db/KEGG/97.0/ko_list
PFAM=/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm

export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/kofamscan/1.3.0/env/bin:/opt/linux/rocky/8.x/x86_64/pkgs/hmmer/3.4/bin:/usr/bin:/bin

BASE=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
WORK=$BASE/work/annot_154
CONF=$WORK/kofam_config.yml
LIST=$BASE/work/cazy_focal/cazy_targets.tsv
PFOUT=$WORK/pfam
KOOUT=$WORK/kofam
TMPD=$WORK/tmp

for B in "$KOFAM" "$HMMSEARCH"; do
  if [ ! -x "$B" ]; then echo "FATAL: not executable: $B"; exit 1; fi
done
for F in "$KOLIST" "$CONF" "$PFAM" "${PFAM}.h3i" "${PFAM}.h3f" "${PFAM}.h3m" "${PFAM}.h3p"; do
  if [ ! -s "$F" ]; then echo "FATAL: missing or empty: $F"; exit 1; fi
done
if [ ! -d "$KOPROF" ]; then echo "FATAL: missing $KOPROF"; exit 1; fi
if [ ! -f "$LIST" ]; then echo "FATAL: missing $LIST"; exit 1; fi

mkdir -p "$PFOUT" "$KOOUT" "$TMPD"

LINE=$(awk -F'\t' -v i="$SLURM_ARRAY_TASK_ID" 'NR>1 && $1==i {print; exit}' "$LIST")
if [ -z "$LINE" ]; then echo "no row for index $SLURM_ARRAY_TASK_ID"; exit 0; fi

GID=$(echo "$LINE" | cut -f2)
FAA=$(echo "$LINE" | cut -f3)
SET=$(echo "$LINE" | cut -f4)

if [ ! -s "$FAA" ]; then echo "FATAL: proteome missing: $FAA"; exit 1; fi
NP=$(grep -c '^>' "$FAA")
echo "index $SLURM_ARRAY_TASK_ID  genome $GID  set $SET  proteins $NP"
if [ "$NP" -lt 300 ]; then echo "FATAL: only $NP proteins"; exit 1; fi

PF=$PFOUT/${GID}.domtbl
if [ -s "$PF" ]; then
  echo "pfam already done: $(grep -vc '^#' "$PF") rows"
else
  echo "=== PFAM 38.2, --cut_ga ==="
  $HMMSEARCH --domtblout "$PF.partial" --cut_ga --cpu 8 -o /dev/null "$PFAM" "$FAA"
  if [ $? -ne 0 ]; then echo "FATAL: pfam failed for $GID"; rm -f "$PF.partial"; exit 1; fi
  mv "$PF.partial" "$PF"
  echo "pfam rows: $(grep -vc '^#' "$PF")"
fi

KO=$KOOUT/${GID}.kofam.tsv
if [ -s "$KO" ]; then
  echo "kofam already done"
else
  echo "=== KOFAMSCAN, per-KO adaptive thresholds ==="
  TD=$TMPD/${GID}
  rm -rf "$TD"; mkdir -p "$TD"
  $KOFAM -c "$CONF" -o "$KO.partial" -p "$KOPROF" -k "$KOLIST" \
    --cpu 8 -f detail-tsv --no-report-unannotated --tmp-dir "$TD" "$FAA"
  RC=$?
  rm -rf "$TD"
  if [ $RC -ne 0 ]; then echo "FATAL: kofamscan exit $RC for $GID"; rm -f "$KO.partial"; exit $RC; fi
  mv "$KO.partial" "$KO"
  echo "kofam rows: $(grep -vc '^#' "$KO")"
  echo "kofam above-threshold: $(awk -F'\t' '$1=="*"' "$KO" | wc -l)"
fi

echo "ANNOT_154_V2_20260806 TASK $SLURM_ARRAY_TASK_ID OK"
