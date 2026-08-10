#!/bin/bash
# The wild SGB representatives are classified against GTDB r220 with GTDB-Tk 2.4.1.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/27_gtdbtk_wild_sgb.sh
# Output: results/gtdbtk_wild_sgb_r220/gtdbtk.bac120.summary.tsv
# GTDBTK_WILD_SGB_V1_20260803
#SBATCH --job-name=gtdbtk_wild
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/gtdbtk_wild_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/gtdbtk_wild_%j.err

# Absolute log paths above. A relative path resolves against the submit
# directory and the job dies with no log (27102139).
# No `set -euo pipefail`: it causes silent failures in SLURM shells here.
# No `module purge`: it breaks the module system in batch shells.

ROOT=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
BATCH=$ROOT/data/tasks/wild_sgb_gtdbtk_batchfile.tsv
OUTDIR=$ROOT/results/gtdbtk_wild_sgb_r220
TMPDIR=$ROOT/work/tmp_gtdbtk_wild

module load gtdbtk/2.4.1

export GTDBTK_DATA_PATH=/srv/projects/db/gtdbtk/220
mkdir -p $OUTDIR
mkdir -p $TMPDIR

echo "host        : $(hostname)"
echo "job         : $SLURM_JOB_ID"
echo "started     : $(date)"
echo "batchfile   : $BATCH"
echo "genomes in  : $(wc -l < $BATCH)"
echo "db          : $GTDBTK_DATA_PATH"
echo "gtdbtk      : $(which gtdbtk)"
gtdbtk --version
echo "-----------------------------------------------------------"

# --pplacer_cpus 8, NEVER 1. At 1 CPU the EHI class-level pass was still
# running at 11h44m and had to be killed (27098673). At 8 on highmem with
# 500 G the same work took 2h53m (27122576).
gtdbtk classify_wf \
  --batchfile $BATCH \
  --out_dir $OUTDIR \
  --cpus 16 \
  --pplacer_cpus 8 \
  --skip_ani_screen \
  --tmpdir $TMPDIR

echo "-----------------------------------------------------------"
echo "gtdbtk exit : $?"
echo "finished    : $(date)"
echo "summary     : $OUTDIR/gtdbtk.bac120.summary.tsv"
if [ -f $OUTDIR/gtdbtk.bac120.summary.tsv ]; then
  echo "rows        : $(( $(wc -l < $OUTDIR/gtdbtk.bac120.summary.tsv) - 1 ))"
fi
echo "GTDBTK_WILD_SGB_FINISHED"
