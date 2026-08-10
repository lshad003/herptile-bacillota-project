#!/bin/bash -l
# Youngblut SGB representatives are classified against GTDB r220.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/27_gtdbtk_youngblut.sh
# Output: results/gtdbtk_youngblut_r220/gtdbtk.bac120.summary.tsv
#SBATCH --job-name=gtdbtk_yb
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/gtdbtk_yb_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/gtdbtk_yb_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --time=24:00:00
#SBATCH --partition=highmem

# SBATCH RATIONALE: classify_wf on 393 Youngblut Clostridia SGB representatives.
# Job 27122576 classified 2,481 genomes in 2h53m at 8 pplacer CPUs and 86 G RSS,
# so 393 should take well under an hour. Same resources because pplacer memory
# scales with the reference tree, not the input count.
#
# --pplacer_cpus 8, NOT 1. At 1 CPU the same work took over 18 hours (27098673).
# r220 to match the herptile MAGs, the EHI arm and the family trees.
#
# NO 'module purge'. NO 'set -ue'.

WORKDIR=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
BATCH=${WORKDIR}/data/youngblut/youngblut_gtdbtk_batchfile.tsv
OUTDIR=${WORKDIR}/results/gtdbtk_youngblut_r220

source /etc/profile.d/modules.sh || true
module load gtdbtk

if ! command -v gtdbtk > /dev/null; then
    echo "[FATAL] gtdbtk not on PATH" >&2
    exit 1
fi
echo "[INFO] gtdbtk: $(command -v gtdbtk)"
gtdbtk --version

export GTDBTK_DATA_PATH=/srv/projects/db/gtdbtk/220
if [ ! -d "${GTDBTK_DATA_PATH}/pplacer" ]; then
    echo "[FATAL] r220 database incomplete" >&2
    exit 1
fi

if [ ! -s "${BATCH}" ]; then
    echo "[FATAL] batchfile missing: ${BATCH}" >&2
    exit 1
fi
N=$(wc -l < ${BATCH})
echo "[INFO] genomes: ${N}"
head -2 ${BATCH}

mkdir -p ${OUTDIR}
export TMPDIR=${WORKDIR}/work/tmp_gtdbtk_yb
mkdir -p ${TMPDIR}

echo "[START] $(date) host=$(hostname)"
gtdbtk classify_wf \
    --batchfile ${BATCH} \
    --out_dir ${OUTDIR} \
    --extension gz \
    --cpus ${SLURM_CPUS_PER_TASK} \
    --pplacer_cpus 8 \
    --skip_ani_screen \
    --tmpdir ${TMPDIR}
RC=$?
echo "[DONE] $(date) rc=${RC}"

if [ ${RC} -ne 0 ]; then
    echo "[FATAL] gtdbtk rc=${RC}" >&2
    exit 1
fi

SUM=${OUTDIR}/gtdbtk.bac120.summary.tsv
if [ -s "${SUM}" ]; then
    echo "[INFO] classified: $(( $(wc -l < ${SUM}) - 1 )) of ${N}"
    echo "[INFO] top families:"
    cut -f2 ${SUM} | tail -n +2 | sed 's/.*;f__//;s/;.*//' | sort | uniq -c | sort -rn | head -8
    echo "[INFO] Ruminococcaceae genera:"
    awk -F'\t' 'NR>1 && $2 ~ /f__Ruminococcaceae/' ${SUM} \
      | cut -f2 | sed 's/.*;g__//;s/;.*//' | sort | uniq -c | sort -rn | head -20
else
    echo "[WARN] no summary at ${SUM}"
fi
echo "GTDBTK_YOUNGBLUT_FINISHED"
