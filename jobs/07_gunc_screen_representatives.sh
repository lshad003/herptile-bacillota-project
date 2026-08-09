#!/bin/bash -l
# Chimerism screening of the SGB representatives against proGenomes 2.1.
#
# Source: ruminococcaceae-agent/jobs/24_gunc_sgb.sh
# Writes: results/gunc_sgb/
#
# GUNC 1.0.6 pins DIAMOND to version 2.0.4 and ships that binary in its own
# module. Loading the diamond module places a later version ahead of it on
# PATH and GUNC aborts.
#
# Two SGBs whose representatives were replaced after exclusion of the
# laboratory-reared arm were screened separately; see
# jobs/09_screen_replacement_representatives.sh.
#SBATCH --job-name=gunc_sgb
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/gunc_sgb_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/gunc_sgb_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=batch

# SBATCH RATIONALE: DIAMOND blastp of 1,171 proteomes (~2.9M proteins) against
# a 13.5 GB reference database. 64G is comfortable, 24 threads is where
# DIAMOND scales. Gene calling is SKIPPED because results/prodigal was called
# with prodigal -p meta, which is GUNC's own default, so --gene_calls saves
# the slowest step.
#
# DO NOT `module load diamond`. GUNC 1.0.6 pins DIAMOND to exactly 2.0.4 and
# ships that binary in its own module bin. Loading the diamond module puts
# 2.1.24 ahead of it on PATH and GUNC aborts:
#   "Diamond version is 2.1.24, not 2.0.4"  (job 27105002)
#
# Log paths are ABSOLUTE. Job 27102139 failed with no logs at all because the
# relative path resolved against the submit directory, which was ~.
#
# NO 'module purge'. NO 'set -ue'.

WORKDIR=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
IN=${WORKDIR}/work/gunc_input
OUT=${WORKDIR}/results/gunc_sgb
DB=/srv/projects/db/GUNC/gunc_db_progenomes2.1.dmnd

source /etc/profile.d/modules.sh || true
module load gunc

if ! command -v gunc > /dev/null; then
    echo "[FATAL] gunc not on PATH" >&2
    exit 1
fi
if ! command -v diamond > /dev/null; then
    echo "[FATAL] diamond not on PATH" >&2
    exit 1
fi
echo "[INFO] gunc:    $(command -v gunc)"
echo "[INFO] diamond: $(command -v diamond)"
gunc --version
DV=$(diamond --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
echo "[INFO] diamond version: ${DV}"
if [ "${DV}" != "2.0.4" ]; then
    echo "[FATAL] GUNC 1.0.6 requires diamond 2.0.4, found ${DV}." >&2
    echo "        Something put another diamond ahead of GUNC's own on PATH." >&2
    exit 1
fi

if [ ! -s "${DB}" ]; then
    echo "[FATAL] GUNC database missing: ${DB}" >&2
    exit 1
fi
echo "[INFO] db: ${DB} ($(du -h ${DB} | cut -f1))"

N=$(ls -1 ${IN}/*.faa 2>/dev/null | wc -l)
if [ "${N}" -lt 1 ]; then
    echo "[FATAL] no .faa in ${IN}" >&2
    exit 1
fi
echo "[INFO] proteomes: ${N}"

mkdir -p ${OUT}
export TMPDIR=${WORKDIR}/work/tmp_gunc
mkdir -p ${TMPDIR}

echo "[START] $(date) host=$(hostname)"
gunc run \
    --input_dir ${IN} \
    --gene_calls \
    --file_suffix .faa \
    --db_file ${DB} \
    --out_dir ${OUT} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --temp_dir ${TMPDIR} \
    --detailed_output
RC=$?
echo "[DONE] $(date) rc=${RC}"

if [ ${RC} -ne 0 ]; then
    echo "[FATAL] gunc rc=${RC}" >&2
    exit 1
fi

SUM=$(ls -1 ${OUT}/GUNC.*.maxCSS_level.tsv 2>/dev/null | head -1)
if [ -s "${SUM}" ]; then
    echo "[INFO] summary: ${SUM}"
    echo "[INFO] rows: $(( $(wc -l < ${SUM}) - 1 )) of ${N}"
    echo "[INFO] pass.GUNC at default CSS 0.45:"
    awk -F'\t' 'NR>1{print $NF}' ${SUM} | sort | uniq -c
    echo "[INFO] taxonomic level of max CSS:"
    awk -F'\t' 'NR>1{print $2}' ${SUM} | sort | uniq -c | sort -rn
else
    echo "[WARN] no maxCSS_level.tsv in ${OUT}"
    ls -1 ${OUT} | head -10
fi
echo "GUNC_FINISHED"
