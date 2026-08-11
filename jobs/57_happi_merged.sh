#!/bin/bash
# The real fit and two label permutations are run as one array.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/39_happi_merged.sh
# Output: work/focal_genus_pangenome/matrices/happi_results_og_bacteria.tsv and permuted controls
# HAPPI_MERGED_V1_20260805
#SBATCH --job-name=happi_og
#SBATCH --partition=stajichlab
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-3
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/happi_og_%A_%a.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/happi_og_%A_%a.err

# Prevalence at ORTHOLOGOUS GROUP level, not raw MMseqs2 cluster level.
# In the unmerged matrix a single gene was split across several clusters,
# so every prevalence was deflated: rpsU sat at 0.76/0.58 and dnaK below
# 1.00 despite being universal, and only 1 cluster reached 95% prevalence
# in the amphibian arm. After merging by eggNOG Bacteria-level OG, dnaK,
# gyrA, sun, apt, prkC and cel all reach 0.98-1.00 in both arms and the
# amphibian core rises from 1 to 350 groups.
#
# Tasks 0-1 are the real comparison at the two merge levels. Tasks 2-3
# repeat the primary on PERMUTED group labels as a null control: if the
# real run returns hits and the permuted runs return none, the method is
# not manufacturing signal.

ROOT=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
case $SLURM_ARRAY_TASK_ID in
  0) TAG=og_bacteria; PERM=0; SEED=0 ;;
  1) TAG=og_narrow;   PERM=0; SEED=0 ;;
  2) TAG=og_bacteria; PERM=1; SEED=101 ;;
  3) TAG=og_bacteria; PERM=1; SEED=202 ;;
esac

echo "host  : $(hostname)"
echo "task  : $SLURM_ARRAY_TASK_ID  tag=$TAG perm=$PERM seed=$SEED"
echo "start : $(date)"
echo "-----------------------------------------------------------"

cd $ROOT
CLUSTERING=$TAG PERMUTE=$PERM PERM_SEED=$SEED Rscript scripts/run_happi_merged.R
RC=$?

echo "-----------------------------------------------------------"
echo "exit  : $RC"
echo "end   : $(date)"
if [ $RC -ne 0 ]; then echo "HAPPI FAILED $TAG perm=$PERM"; exit $RC; fi
echo "HAPPI_MERGED_FINISHED $TAG perm=$PERM"
