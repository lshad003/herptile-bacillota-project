#!/bin/bash
# EHI newt genomes are classified against GTDB r220.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/31_gtdbtk_ehi_amphibian.sh
# Output: results/gtdbtk_ehi_amphibian_r220/gtdbtk.bac120.summary.tsv
# GTDBTK_EHI_AMPH_V1_20260804
#SBATCH --job-name=gtdbtk_amph
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/gtdbtk_amph_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/gtdbtk_amph_%j.err

# Identical settings to job 27138736 (718 wild SGB reps, 1h21m) so this arm
# is classified against the same r220 database with the same parameters.
# --pplacer_cpus 8, NEVER 1: at 1 CPU the EHI class-level pass ran 11h44m
# and had to be killed (27098673).
# Absolute log paths above; a relative path resolves against the submit
# directory and the job dies with no log (27102139).

ROOT=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
BATCH=$ROOT/data/tasks/ehi_amphibian_gtdbtk_batchfile.tsv
OUTDIR=$ROOT/results/gtdbtk_ehi_amphibian_r220
TMPDIR=$ROOT/work/tmp_gtdbtk_amph

module load gtdbtk/2.4.1
export GTDBTK_DATA_PATH=/srv/projects/db/gtdbtk/220
mkdir -p $OUTDIR
mkdir -p $TMPDIR

echo "host       : $(hostname)"
echo "job        : $SLURM_JOB_ID"
echo "started    : $(date)"
echo "batchfile  : $BATCH"
echo "genomes in : $(wc -l < $BATCH)"
echo "db         : $GTDBTK_DATA_PATH"
gtdbtk --version
echo "-----------------------------------------------------------"

gtdbtk classify_wf \
  --batchfile $BATCH \
  --out_dir $OUTDIR \
  --cpus 16 \
  --pplacer_cpus 8 \
  --skip_ani_screen \
  --tmpdir $TMPDIR

RC=$?
echo "-----------------------------------------------------------"
echo "gtdbtk exit : $RC"
if [ $RC -ne 0 ]; then echo "GTDBTK FAILED"; exit $RC; fi

SUM=$OUTDIR/gtdbtk.bac120.summary.tsv
if [ ! -s $SUM ]; then echo "NO SUMMARY PRODUCED"; exit 1; fi

echo "finished   : $(date)"
echo "summary    : $SUM"
echo "rows       : $(( $(wc -l < $SUM) - 1 ))"
echo "Ruminococcaceae: $(grep -c 'f__Ruminococcaceae' $SUM)"
echo "GTDBTK_EHI_AMPH_FINISHED"
