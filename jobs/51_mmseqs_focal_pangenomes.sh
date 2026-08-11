#!/bin/bash -l
# Focal proteomes are clustered with MMseqs2.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/21_mmseqs_focal_pangenomes.sh
# Output: work/focal_genus_pangenome/mmseqs/
#SBATCH --job-name=mmseqs_focal
#SBATCH --output=logs/mmseqs_focal_%j.out
#SBATCH --error=logs/mmseqs_focal_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=batch

# SBATCH RATIONALE: three easy-cluster runs, 110k to 161k proteins each.
# Job 26995090 did 1,127,761 proteins in 9m58s using 2.11 G at 16 threads, so
# each of these is roughly 14% of that. Small enough to be login-node safe in
# principle, but MMseqs2 on this project has killed a login node before and
# 16G/1h costs nothing.
#
# NO 'module purge' (broke the module system in batch shells on Jul 30).
# NO 'set -ue' (silent SLURM failures). Errors reported explicitly.
# Parameters are covmode1, chosen by the rule fixed before the numbers were
# seen: most high-prevalence clusters (88) at 1.106 proteins per genome.
# Each genus writes into its own directory. Nothing touches work/pangenome.

source /etc/profile.d/modules.sh || true
module load mmseqs2/15-6f452

BASE=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/work/focal_pangenomes

if ! command -v mmseqs > /dev/null; then
    echo "[FATAL] mmseqs not on PATH" >&2
    exit 1
fi
echo "[INFO] using $(command -v mmseqs)"

run_one () {
    GEN=$1
    D="${BASE}/${GEN}"
    FAA="${D}/${GEN}_proteins.faa"

    if [ ! -s "${FAA}" ]; then
        echo "[FATAL] missing input ${FAA}" >&2
        return 1
    fi
    N=$(grep -c '>' "${FAA}")
    echo "[INFO] ${GEN} input proteins: ${N}"

    DUP=$(grep '>' "${FAA}" | sort | uniq -d | head -3)
    if [ -n "${DUP}" ]; then
        echo "[FATAL] ${GEN} duplicate protein ids:" >&2
        echo "${DUP}" >&2
        return 1
    fi

    cd "${D}" || return 1
    rm -rf tmp clu_all_seqs.fasta clu_cluster.tsv clu_rep_seq.fasta
    mkdir -p tmp

    echo "[START] ${GEN} $(date)"
    mmseqs easy-cluster "${GEN}_proteins.faa" clu tmp \
        --min-seq-id 0.5 -c 0.8 --cov-mode 1 \
        --threads ${SLURM_CPUS_PER_TASK}
    RC=$?
    rm -rf tmp

    if [ ${RC} -ne 0 ]; then
        echo "[FATAL] ${GEN} mmseqs rc=${RC}" >&2
        return 1
    fi
    if [ ! -s clu_cluster.tsv ]; then
        echo "[FATAL] ${GEN} no clu_cluster.tsv" >&2
        return 1
    fi

    LINES=$(wc -l < clu_cluster.tsv)
    CLUS=$(cut -f1 clu_cluster.tsv | sort -u | wc -l)
    echo "[DONE] ${GEN} $(date)"
    echo "[${GEN}] lines: ${LINES} (expected ${N})"
    echo "[${GEN}] clusters: ${CLUS}"
    if [ "${LINES}" != "${N}" ]; then
        echo "[WARN] ${GEN} line count does not match input protein count" >&2
    fi
    return 0
}

FAIL=0
for G in UBA866 Anaerotruncus Angelakisella; do
    run_one "${G}" || FAIL=1
done

if [ ${FAIL} -ne 0 ]; then
    echo "[FATAL] at least one genus failed" >&2
    exit 1
fi
echo "ALL_FOCAL_PANGENOMES_FINISHED"
