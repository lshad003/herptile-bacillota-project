#!/bin/bash -l
# barrnap and tRNAscan-SE on 1,247 GTDB r220 Ruminococcaceae references.
# Same settings as the catalog arm, so recovery rates are comparable.
# Source: ruminococcaceae-agent/jobs/51_rrna_trna_refs.sh
# Output: results/rrna_trna_refs/{barrnap,trnascan}/

# RRNA_TRNA_REFS_V1_20260806
#SBATCH --job-name=rt_refs
#SBATCH --partition=intel
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=0-62%30
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/rt_refs_%A_%a.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/rt_refs_%A_%a.err

# COMPARISON ARM for the MIMAG rRNA+tRNA screen. On the 1,171 SGB
# representatives only 3 (0.4%) carry 5S+16S+23S plus tRNAs for >=18 amino
# acids, and ZERO Ruminococcaceae do. Without a reference arm that reads as
# a defect specific to this catalog. 91.9% of GTDB r220 Ruminococcaceae are
# themselves MAGs (R7), so if they score similarly the correct statement is
# that MIMAG high quality is not attainable for this family from short-read
# metagenome assembly, which is a claim about the method not the data.
#
# CHATINDEX R5 warns that an earlier "70-85% for published human gut MAG
# catalogs" comparison was UNSOURCED. This measures it instead.
#
# Genomes are the 1,247 already staged for the reference antiSMASH run
# (jobs/46_antismash_refs.sh), gzipped, one per line in bgc_ref_input_list.txt.
#
# barrnap/0.9 and trnascan-se modules load miniconda3, which fails in a
# batch shell. Both are called by absolute path and PATH is set explicitly.
# 1,247 genomes at CHUNK=20 needs 63 tasks (0-62). Task 62 handles 1241-1247.

ROOT=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
LIST=$ROOT/work/bgc_refs/bgc_ref_input_list.txt
OUT=$ROOT/results/rrna_trna_refs
SCRATCH=$ROOT/work/rt_refs_scratch
CHUNK=20

BARRNAP=/opt/linux/rocky/8.x/x86_64/pkgs/barrnap/0.9/bin/barrnap
TRNASCAN=/opt/linux/rocky/8.x/x86_64/pkgs/trnascan-se/2.0.12/bin/tRNAscan-SE
export PERL5LIB=/opt/linux/rocky/8.x/x86_64/pkgs/trnascan-se/2.0.12/bin:/opt/linux/rocky/8.x/x86_64/pkgs/trnascan-se/2.0.12/lib/tRNAscan-SE:$PERL5LIB
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/barrnap/0.9/bin:/opt/linux/rocky/8.x/x86_64/pkgs/trnascan-se/2.0.12/bin:$PATH

mkdir -p $OUT/barrnap $OUT/trnascan $SCRATCH

echo "host  : $(hostname)"
echo "task  : $SLURM_ARRAY_TASK_ID"
echo "start : $(date)"
if [ ! -x "$BARRNAP" ];  then echo "BARRNAP NOT EXECUTABLE: $BARRNAP"; exit 1; fi
if [ ! -x "$TRNASCAN" ]; then echo "TRNASCAN NOT EXECUTABLE: $TRNASCAN"; exit 1; fi
if [ ! -s "$LIST" ];     then echo "MISSING LIST: $LIST"; exit 1; fi
$BARRNAP --version 2>&1 | head -1
$TRNASCAN -h 2>&1 | head -2

TOTAL=$(wc -l < $LIST)
START=$(( SLURM_ARRAY_TASK_ID * CHUNK + 1 ))
END=$(( START + CHUNK - 1 ))
if [ $START -gt $TOTAL ]; then echo "beyond $TOTAL genomes, nothing to do"; exit 0; fi
echo "genomes total: $TOTAL, this task handles $START-$END"
echo "-----------------------------------------------------------"

N=0; NB=0; NT=0
for i in $(seq $START $END); do
  [ $i -gt $TOTAL ] && break
  GZ=$(sed -n "${i}p" $LIST)
  [ -z "$GZ" ] && continue
  [ -s "$GZ" ] || { echo "MISSING INPUT: $GZ"; continue; }
  B=$(basename $GZ); B=${B%.fna.gz}; B=${B%.fna}

  if [ -s $OUT/barrnap/${B}.gff ] && [ -s $OUT/trnascan/${B}.txt ]; then
    NB=$((NB+1)); NT=$((NT+1)); N=$((N+1)); continue
  fi

  FA=$SCRATCH/${B}.${SLURM_ARRAY_JOB_ID}.${SLURM_ARRAY_TASK_ID}.fna
  zcat $GZ > $FA
  if [ ! -s "$FA" ]; then echo "EMPTY AFTER ZCAT: $B"; rm -f $FA; continue; fi

  if [ ! -s $OUT/barrnap/${B}.gff ]; then
    $BARRNAP --kingdom bac --threads 4 --quiet $FA > $OUT/barrnap/${B}.gff 2>/dev/null
  fi
  [ -s $OUT/barrnap/${B}.gff ] && NB=$((NB+1))

  if [ ! -s $OUT/trnascan/${B}.txt ]; then
    $TRNASCAN -B -q -o $OUT/trnascan/${B}.txt $FA > /dev/null 2>&1
  fi
  [ -s $OUT/trnascan/${B}.txt ] && NT=$((NT+1))

  rm -f $FA
  N=$((N+1))
done

echo "genomes attempted: $N, barrnap outputs: $NB, trnascan outputs: $NT"
if [ $NB -eq 0 ] || [ $NT -eq 0 ]; then
  echo "ONE OR BOTH TOOLS PRODUCED NOTHING"
  exit 1
fi
echo "end   : $(date)"
echo "RRNA_TRNA_REFS_FINISHED task $SLURM_ARRAY_TASK_ID"
