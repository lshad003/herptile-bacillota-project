#!/bin/bash -l
# Detection of rRNA genes and tRNAs in the 1,171 SGB representatives.
#
# Source: ruminococcaceae-agent/jobs/45_rrna_trna.sh
# Note:   a resume job, jobs/48_rrna_trna_resume.sh, was run afterwards to
#         cover tasks that did not complete in the first submission. Both
#         write into the same output directories.
# Writes: results/rrna_trna/barrnap/, results/rrna_trna/trnascan/
#
# barrnap writes a GFF containing only a version header when nothing is
# found. Such a file is non-empty, so presence is determined by counting
# feature lines by rRNA type rather than by testing whether the file exists.
# RRNA_TRNA_V2_20260805
#SBATCH --job-name=rrna_trna
#SBATCH --partition=short
#SBATCH --time=01:45:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=0-58%30
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/rrna_trna_%A_%a.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/rrna_trna_%A_%a.err

# V1 (27245223, 27245225) failed on every task in 2-4 seconds: the
# barrnap/0.9 module does `module load miniconda3`, which fails in a batch
# shell the same way dRep did in job 27148950. Both binaries are now called
# by absolute path. trnascan-se only prepends PATH and would have worked
# either way, but is resolved too for consistency.
#
# CHATINDEX R1 forbids "MIMAG high quality" because MIMAG also requires 5S,
# 16S and 23S rRNA plus tRNAs for at least 18 of 20 amino acids, never
# checked. This covers the 1,171 SGB representatives. The 49.3%
# near-complete figure in Results 3.1 is per-MAG across all 2,229 and would
# need a second run on that set to use MIMAG language.

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
