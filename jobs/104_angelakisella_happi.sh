#!/bin/bash
# Real fit and two label permutations.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/angelakisella_happi.sh
# Output: results/happi_angelakisella_og.tsv, results/happi_angelakisella_og_perm101.tsv, results/happi_angelakisella_og_perm202.tsv
#SBATCH --job-name=happi_angel
#SBATCH --partition=stajichlab
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-2
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/happi_angel_%A_%a.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/happi_angel_%A_%a.err

# happi on the Angelakisella OG matrix, 31 amphibian vs 26 non-amphibian,
# completeness-aware. Task 0 = real labels; tasks 1-2 = permuted labels
# (seeds 101, 202, the project's standard permutation seeds): if the
# real run returns hits and both permutations return none, the method is
# not manufacturing signal. Parameters identical to jobs/39_happi_merged.sh.

ROOT=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
case $SLURM_ARRAY_TASK_ID in
  0) PERM=0; SEED=0 ;;
  1) PERM=1; SEED=101 ;;
  2) PERM=1; SEED=202 ;;
esac
echo "host  : $(hostname)"
echo "task  : $SLURM_ARRAY_TASK_ID  perm=$PERM seed=$SEED"
echo "start : $(date)"
echo "-----------------------------------------------------------"
cd $ROOT
PERMUTE=$PERM PERM_SEED=$SEED Rscript scripts/run_happi_angelakisella.R
RC=$?
echo "-----------------------------------------------------------"
echo "exit  : $RC"
echo "end   : $(date)"
if [ $RC -ne 0 ]; then echo "HAPPI FAILED perm=$PERM seed=$SEED"; exit $RC; fi
echo "HAPPI_ANGELAKISELLA_JOB_FINISHED perm=$PERM seed=$SEED"
# ANGELAKISELLA_HAPPI_JOB_V1
