#!/bin/bash -l
# The classify step is rerun against the completed alignment, producing the EHI assignments used here.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/26_gtdbtk_classify_resume.sh
# Output: results/gtdbtk_ehi_r220_classify/gtdbtk.bac120.summary.tsv
#SBATCH --job-name=gtdbtk_res
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/gtdbtk_resume_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/gtdbtk_resume_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --partition=highmem

# SBATCH RATIONALE: resumes ONLY the classify step of job 27098673, reusing the
# completed identify/ (16 G) and align/ (240 M) output. Those two stages took
# ~85 minutes; the expensive part is pplacer.
#
# WHY THIS DIFFERS FROM 27098673:
#   --pplacer_cpus 1 -> 8. The first pplacer pass took 4h54m at 1 CPU and the
#   second was still running at 11h when the 24h wall limit approached.
#   pplacer memory scales with pplacer_cpus, hence 500G rather than 160G.
#   Time raised to 48h; highmem allows 30 days, so the wall limit is no longer
#   the binding constraint.
#
# highmem nodes: 7 x 922 GB, 32+ CPUs. 500G schedules comfortably.
#
# DOES NOT TOUCH results/gtdbtk_ehi_r220. Writes to a new directory and reads
# the old one through --align_dir, so identify/ and align/ cannot be lost.
#
# NO 'module purge'. NO 'set -ue'.

WORKDIR=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
OLD=${WORKDIR}/results/gtdbtk_ehi_r220
NEW=${WORKDIR}/results/gtdbtk_ehi_r220_classify
BATCH=${WORKDIR}/data/tasks/ehi_gtdbtk_batchfile.tsv

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

echo "[INFO] classify flags available in this version:"
gtdbtk classify --help 2>&1 | grep -E 'align_dir|batchfile|genome_dir|skip_ani_screen|pplacer_cpus|mash_db' | head -12

if [ ! -d "${OLD}/align" ]; then
    echo "[FATAL] no align dir at ${OLD}/align" >&2
    exit 1
fi
echo "[INFO] contents of ${OLD}/align:"
ls -la ${OLD}/align | head -12

MSA=$(ls -1 ${OLD}/align/*user_msa* 2>/dev/null | head -1)
if [ -z "${MSA}" ]; then
    echo "[FATAL] no user_msa file found in ${OLD}/align" >&2
    echo "        align may not have completed. Do not resume." >&2
    exit 1
fi
echo "[INFO] user MSA: ${MSA} ($(du -h ${MSA} | cut -f1))"

if [ ! -s "${BATCH}" ]; then
    echo "[FATAL] batchfile missing: ${BATCH}" >&2
    exit 1
fi
echo "[INFO] genomes in batchfile: $(wc -l < ${BATCH})"

mkdir -p ${NEW}
export TMPDIR=${WORKDIR}/work/tmp_gtdbtk_resume
mkdir -p ${TMPDIR}

echo "[START] $(date) host=$(hostname)"
gtdbtk classify \
    --batchfile ${BATCH} \
    --align_dir ${OLD} \
    --out_dir ${NEW} \
    --extension gz \
    --cpus ${SLURM_CPUS_PER_TASK} \
    --pplacer_cpus 8 \
    --skip_ani_screen \
    --tmpdir ${TMPDIR}
RC=$?
echo "[DONE] $(date) rc=${RC}"

if [ ${RC} -ne 0 ]; then
    echo "[FATAL] gtdbtk classify rc=${RC}" >&2
    exit 1
fi

SUM=${NEW}/gtdbtk.bac120.summary.tsv
if [ -s "${SUM}" ]; then
    echo "[INFO] classified: $(( $(wc -l < ${SUM}) - 1 ))"
    echo "[INFO] top families:"
    cut -f2 ${SUM} | tail -n +2 | sed 's/.*;f__/f__/;s/;.*//' | sort | uniq -c | sort -rn | head -8
    echo "[INFO] focal genera:"
    for G in Anaerotruncus UBA866 Paludihabitans Angelakisella Ruthenibacterium Fimivivens; do
        C=$(cut -f2 ${SUM} | grep -c ";g__${G};")
        echo "    ${G}: ${C}"
    done
else
    echo "[WARN] no summary at ${SUM}"
    ls -1 ${NEW} | head
fi
echo "GTDBTK_RESUME_FINISHED"
