#!/bin/bash -l
# EHI genomes are classified against GTDB r220. Identify and align completed; classify did not.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/23_gtdbtk_ehi_r220.sh
# Output: results/gtdbtk_ehi_r220/identify/ and align/
#SBATCH --job-name=gtdbtk_ehi
#SBATCH --output=logs/gtdbtk_ehi_%j.out
#SBATCH --error=logs/gtdbtk_ehi_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=24:00:00
#SBATCH --partition=highmem

# SBATCH RATIONALE: GTDB-Tk classify_wf on 2,481 EHI genomes. pplacer is the
# memory bottleneck and scales with the reference tree, not the input count.
# highmem nodes have 900+ GB so 160G is comfortable. --pplacer_cpus 1 keeps it
# to one pplacer process; pplacer memory scales with pplacer_cpus, so if this
# OOMs raise memory, never threads.
#
# r220 NOT r226: GTDB-Tk 2.4.1 ships for r220, and r226 needs 2.5+. r220 also
# matches the herptile MAGs, the reference set, the family trees and the
# candidate-genus table, so both arms of the comparison stay on one release.
#
# NO 'module purge'. NO 'set -ue'.

WORKDIR=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
BATCH=${WORKDIR}/data/tasks/ehi_gtdbtk_batchfile.tsv
OUTDIR=${WORKDIR}/results/gtdbtk_ehi_r220

source /etc/profile.d/modules.sh || true
module load gtdbtk

if ! command -v gtdbtk > /dev/null; then
    echo "[FATAL] gtdbtk not on PATH" >&2
    exit 1
fi
echo "[INFO] gtdbtk: $(command -v gtdbtk)"
gtdbtk --version

export GTDBTK_DATA_PATH=/srv/projects/db/gtdbtk/220
echo "[INFO] GTDBTK_DATA_PATH=${GTDBTK_DATA_PATH}"
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

mkdir -p ${OUTDIR}
export TMPDIR=${WORKDIR}/work/tmp_gtdbtk
mkdir -p ${TMPDIR}

echo "[START] $(date) host=$(hostname)"
gtdbtk classify_wf \
    --batchfile ${BATCH} \
    --out_dir ${OUTDIR} \
    --extension gz \
    --cpus ${SLURM_CPUS_PER_TASK} \
    --pplacer_cpus 1 \
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
    cut -f2 ${SUM} | tail -n +2 | sed 's/.*;f__/f__/;s/;.*//' | sort | uniq -c | sort -rn | head -8
    echo "[INFO] focal genera:"
    for G in Anaerotruncus UBA866 Angelakisella Ruthenibacterium Fimivivens; do
        C=$(cut -f2 ${SUM} | grep -c ";g__${G};")
        echo "    ${G}: ${C}"
    done
else
    echo "[WARN] no summary at ${SUM}"
fi
echo "GTDBTK_EHI_FINISHED"
