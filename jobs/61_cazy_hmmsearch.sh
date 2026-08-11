#!/bin/bash
# Proteomes are searched against dbCAN v13.0 with HMMER.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/63_cazy_hmmsearch.sh
# Output: work/cazy_focal/domtbl/
#SBATCH --job-name=cazy154
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --array=0-153%40
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/cazy_%A_%a.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/cazy_%A_%a.err

HMMSEARCH=/opt/linux/rocky/8.x/x86_64/pkgs/hmmer/3.4/bin/hmmsearch
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/hmmer/3.4/bin:/usr/bin:/bin

BASE=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
WORK=$BASE/work/cazy_focal
LIST=$WORK/cazy_targets.tsv
HMMDB=/srv/projects/db/CAZY/CAZyDB/v13.0/dbCAN-HMMdb-V13.txt
OUTD=$WORK/domtbl

if [ ! -x "$HMMSEARCH" ]; then
  echo "FATAL: hmmsearch not executable at $HMMSEARCH"
  exit 1
fi
for EXT in h3f h3i h3m h3p; do
  if [ ! -s "${HMMDB}.${EXT}" ]; then
    echo "FATAL: pressed index ${HMMDB}.${EXT} missing or empty"
    exit 1
  fi
done
if [ ! -f "$LIST" ]; then
  echo "FATAL: target list missing: $LIST"
  exit 1
fi

mkdir -p "$OUTD"

LINE=$(awk -F'\t' -v i="$SLURM_ARRAY_TASK_ID" 'NR>1 && $1==i {print; exit}' "$LIST")
if [ -z "$LINE" ]; then
  echo "no row for index $SLURM_ARRAY_TASK_ID"
  exit 0
fi

GID=$(echo "$LINE" | cut -f2)
FAA=$(echo "$LINE" | cut -f3)
SET=$(echo "$LINE" | cut -f4)
OUT=$OUTD/${GID}.domtbl

if [ -s "$OUT" ]; then
  echo "already done: $OUT"
  exit 0
fi
if [ ! -s "$FAA" ]; then
  echo "FATAL: proteome missing or empty: $FAA"
  exit 1
fi

NP=$(grep -c '^>' "$FAA")
echo "index $SLURM_ARRAY_TASK_ID  genome $GID  set $SET  proteins $NP"
if [ "$NP" -lt 300 ]; then
  echo "FATAL: only $NP proteins in $FAA"
  exit 1
fi

$HMMSEARCH --domtblout "$OUT.partial" -E 1e-10 --cpu 4 -o /dev/null "$HMMDB" "$FAA"
RC=$?
if [ $RC -ne 0 ]; then
  echo "FATAL: hmmsearch exit $RC for $GID"
  rm -f "$OUT.partial"
  exit $RC
fi

mv "$OUT.partial" "$OUT"
NHIT=$(grep -vc '^#' "$OUT")
echo "wrote $OUT with $NHIT domain rows"
echo "CAZY_HMMSEARCH_V1_20260806 TASK $SLURM_ARRAY_TASK_ID OK"
