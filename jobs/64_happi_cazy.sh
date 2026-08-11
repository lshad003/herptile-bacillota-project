#!/bin/bash
# The real fit and two label permutations are run as one array.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/65_happi_cazy_v2.sh
# Output: results/happi_cazy_focal.tsv, happi_cazy_focal_perm101.tsv, perm202.tsv
#SBATCH --job-name=happicazy
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/happicazy2_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/happicazy2_%j.err

RSCRIPT=/opt/linux/rocky/8.x/x86_64/pkgs/R/4.5.2/bin/Rscript
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/R/4.5.2/bin:/usr/bin:/bin
export R_LIBS_USER=/rhome/lshad003/R/x86_64-pc-linux-gnu-library/4.5

SCRIPT=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/run_happi_cazy_v2.R

if [ ! -x "$RSCRIPT" ]; then
  echo "FATAL: Rscript not executable at $RSCRIPT"
  exit 1
fi
if [ ! -d "$R_LIBS_USER/happi" ]; then
  echo "FATAL: happi not found in $R_LIBS_USER"
  exit 1
fi
if ! grep -q RUN_HAPPI_CAZY_V2_20260806_END "$SCRIPT"; then
  echo "FATAL: end marker missing from $SCRIPT"
  exit 1
fi

echo "script: $SCRIPT  PERMUTE=${PERMUTE:-0} PERM_SEED=${PERM_SEED:-0}"
$RSCRIPT "$SCRIPT"
