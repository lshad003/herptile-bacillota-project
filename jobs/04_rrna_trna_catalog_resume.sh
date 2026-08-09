#!/bin/bash -l
# Resume job for rRNA and tRNA detection in the SGB representatives.
#
# Source: ruminococcaceae-agent/jobs/48_rrna_trna_resume.sh
# Writes: results/rrna_trna/barrnap/, results/rrna_trna/trnascan/
#
# Tasks that did not complete under jobs/03_rrna_trna_catalog.sh were rerun
# here. Output directories are shared between the two jobs, so completeness
# was confirmed by counting output files rather than by job exit state.
# RRNA_TRNA_V3_20260805
#SBATCH --job-name=rrna_trna
#SBATCH --partition=intel
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=1-58%30
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/rrna_trna_%A_%a.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/rrna_trna_%A_%a.err

# V2 (27245543) ran task 0 ONLY, covering genomes 1-20 of 1,171. It
# completed cleanly in 01:12:04, so the logic and the chunk size are sound;
# tasks 1-58 were never submitted. This resumes them.
#
# CHUNK stays at 20 so the task-index-to-genome mapping matches the
# completed task 0. The per-genome guards make a rerun of any finished
# genome a no-op, but task 0 is excluded anyway.
#
# Partition moved off short: 01:12:04 against a 01:45:00 cap left no room
# for a larger-than-average genome.
#
# barrnap/0.9 and trnascan-se modules load miniconda3, which fails in a
# batch shell. Both are called by absolute path and PATH is set explicitly
# so nothing depends on the submitting shell.

ROOT=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
GDIR=$ROOT/results/drep_herptile_95ani_2229/dereplicated_genomes
OUT=$ROOT/results/rrna_trna
CHUNK=20

BARRNAP=/opt/linux/rocky/8.x/x86_64/pkgs/barrnap/0.9/bin/barrnap
TRNASCAN=/opt/linux/rocky/8.x/x86_64/pkgs/trnascan-se/2.0.12/bin/tRNAscan-SE
export PERL5LIB=/opt/linux/rocky/8.x/x86_64/pkgs/trnascan-se/2.0.12/bin:/opt/linux/rocky/8.x/x86_64/pkgs/trnascan-se/2.0.12/lib/tRNAscan-SE:$PERL5LIB
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/barrnap/0.9/bin:/opt/linux/rocky/8.x/x86_64/pkgs/trnascan-se/2.0.12/bin:$PATH

mkdir -p $OUT/barrnap $OUT/trnascan

echo "host  : $(hostname)"
echo "task  : $SLURM_ARRAY_TASK_ID"
echo "start : $(date)"
if [ ! -x "$BARRNAP" ];  then echo "BARRNAP NOT EXECUTABLE: $BARRNAP"; exit 1; fi
if [ ! -x "$TRNASCAN" ]; then echo "TRNASCAN NOT EXECUTABLE: $TRNASCAN"; exit 1; fi
$BARRNAP --version 2>&1 | head -1
$TRNASCAN -h 2>&1 | head -2

LIST=/tmp/rt_list_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt
ls $GDIR/*.fa > $LIST 2>/dev/null
[ -s $LIST ] || ls $GDIR/*.fna > $LIST 2>/dev/null
TOTAL=$(wc -l < $LIST)
START=$(( SLURM_ARRAY_TASK_ID * CHUNK + 1 ))
END=$(( START + CHUNK - 1 ))
if [ $START -gt $TOTAL ]; then echo "beyond $TOTAL genomes, nothing to do"; exit 0; fi
echo "genomes total: $TOTAL, this task handles $START-$END"
echo "-----------------------------------------------------------"

N=0; NB=0; NT=0
for i in $(seq $START $END); do
  [ $i -gt $TOTAL ] && break
  F=$(sed -n "${i}p" $LIST)
  [ -z "$F" ] && continue
  B=$(basename $F); B=${B%.fa}; B=${B%.fna}

  if [ ! -s $OUT/barrnap/${B}.gff ]; then
    $BARRNAP --kingdom bac --threads 4 --quiet $F > $OUT/barrnap/${B}.gff 2>/dev/null
  fi
  [ -s $OUT/barrnap/${B}.gff ] && NB=$((NB+1))

  if [ ! -s $OUT/trnascan/${B}.txt ]; then
    $TRNASCAN -B -q -o $OUT/trnascan/${B}.txt $F > /dev/null 2>&1
  fi
  [ -s $OUT/trnascan/${B}.txt ] && NT=$((NT+1))
  N=$((N+1))
done
rm -f $LIST

echo "genomes attempted: $N, barrnap outputs: $NB, trnascan outputs: $NT"
if [ $NB -eq 0 ] || [ $NT -eq 0 ]; then
  echo "ONE OR BOTH TOOLS PRODUCED NOTHING"
  exit 1
fi
echo "end   : $(date)"
echo "RRNA_TRNA_FINISHED task $SLURM_ARRAY_TASK_ID"
# SENTINEL_END
