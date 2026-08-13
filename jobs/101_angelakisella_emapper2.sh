#!/bin/bash
# eggNOG annotation of the two additional genomes.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/angelakisella_emapper_gap2.sh
# Output: results/eggnog/ and results/eggnog_refs/ annotation files
#SBATCH --job-name=angel_emap2
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/angel_emap2_%A_%a.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/angel_emap2_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --partition=batch
#SBATCH --array=1-2%2

# Annotates the 2 extra genomes for the 61-tip Angelakisella
# neighborhood (Youngblut member + 1 Heteroruminococcus reference).
# Environment per the documented eggnog-mapper batch-shell fix:
# absolute binary path, explicit EGGNOG_DATA_DIR, no module load.
# Parameters match jobs/09b_eggnog_resume.sh.

set -u

BASE=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
STAGE=$BASE/work/angelakisella_emapper_gap2
TASKS=$STAGE/tasks.tsv
EMAPPER=/opt/linux/rocky/8.x/x86_64/pkgs/eggnog-mapper/2.1.9/env/bin/emapper.py
export EGGNOG_DATA_DIR=/srv/projects/db/eggNOG/LATEST
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/eggnog-mapper/2.1.9/env/bin:$PATH
DMND_DB=$EGGNOG_DATA_DIR/eggnog_proteins.dmnd

if [ ! -f "$EMAPPER" ]; then echo "STOP: emapper not found: $EMAPPER" >&2; exit 1; fi
if [ ! -d "$EGGNOG_DATA_DIR" ]; then echo "STOP: db dir not found" >&2; exit 1; fi
if [ ! -f "$DMND_DB" ]; then echo "STOP: dmnd db not found: $DMND_DB" >&2; exit 1; fi

NTASKS=$(($(wc -l < "$TASKS") - 1))
if [[ "$NTASKS" -ne 2 ]]; then
    echo "STOP: task list has $NTASKS tasks, array sized for 2" >&2
    exit 1
fi

LINE=$(awk -F'\t' -v i="$SLURM_ARRAY_TASK_ID" '$1==i {print $2"\t"$3}' "$TASKS")
GID=$(echo "$LINE" | cut -f1)
OUTSUB=$(echo "$LINE" | cut -f2)

if [[ -z "$GID" || -z "$OUTSUB" ]]; then
    echo "STOP: no task for index $SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi

FASTA=$STAGE/faa/$GID.faa
OUTDIR=$BASE/results/$OUTSUB
OUT_ANNO=$OUTDIR/$GID.emapper.annotations

if [[ -s "$OUT_ANNO" ]]; then
    echo "[INFO] already done: $GID"
    exit 0
fi
if [[ ! -s "$FASTA" ]]; then
    echo "[ERROR] missing: $FASTA" >&2
    exit 1
fi

export TMPDIR="/tmp/${USER}.angelemap2.${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"

echo "host: $(hostname)  task: $SLURM_ARRAY_TASK_ID  genome: $GID"
$EMAPPER --version

$EMAPPER \
    -i "$FASTA" --itype proteins \
    -o "$GID" --output_dir "$OUTDIR" \
    --cpu ${SLURM_CPUS_PER_TASK} \
    -m diamond \
    --data_dir "$EGGNOG_DATA_DIR" \
    --dmnd_db "$DMND_DB" \
    --no_file_comments \
    --report_no_hits \
    --sensmode sensitive \
    --evalue 1e-3 \
    --temp_dir "$TMPDIR" \
    --scratch_dir "$TMPDIR" \
    --override
RC=$?

rm -rf "$TMPDIR"
if [ $RC -ne 0 ]; then echo "[ERROR] emapper exit $RC for $GID" >&2; exit $RC; fi
if [ ! -s "$OUT_ANNO" ]; then echo "[ERROR] no annotations file for $GID" >&2; exit 1; fi
echo "[DONE] $(date) $GID  rows: $(grep -vc '^#' "$OUT_ANNO")"
# ANGELAKISELLA_EMAPPER_GAP2_JOB_V1
