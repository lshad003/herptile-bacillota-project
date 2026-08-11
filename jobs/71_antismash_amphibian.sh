#!/bin/bash
# Amphibian genomes are searched with antiSMASH. The metagenomic gene finder is used, since the assemblies are fragmented.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/41_antismash_array.sh
# Output: work/bgc/antismash/
# ANTISMASH_ARRAY_V2_20260805
#SBATCH --job-name=asmash
#SBATCH --partition=stajichlab
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-999%60
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/asmash_%A_%a.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/asmash_%A_%a.err

# V1 (27241732, 27241759) failed on every task with exit 2:
#   argument --genefinding-tool: invalid choice: 'prodigal-meta'
# antiSMASH 7 accepts glimmerhmm, prodigal, prodigal-m, none, error.
# prodigal-m is the metagenomic mode, which is what fragmented MAGs need.
#
# 1,155 genomes: 718 wild SGB representatives and 437 EHI newt Bacillota_A.
# Array max is 1000, so this is submitted twice, OFFSET=0 for tasks 0-999
# then OFFSET=1000 with --array=0-154.
#
# The antismash module activates a conda env and `module load` in a batch
# shell gave CondaError with exit 127 in job 27148950, so the binary is
# called by absolute path.

ROOT=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
WORK=$ROOT/work/bgc
LIST=$WORK/bgc_input_list.txt
OUTBASE=$WORK/antismash
SCRATCH=$WORK/scratch
ASMASH=/opt/linux/rocky/8.x/x86_64/pkgs/antismash/7.1.0/bin/antismash

OFFSET=${OFFSET:-0}
IDX=$(( SLURM_ARRAY_TASK_ID + OFFSET ))
LINE=$(( IDX + 1 ))
TOTAL=$(wc -l < $LIST)

if [ $IDX -ge $TOTAL ]; then
  echo "index $IDX beyond input list ($TOTAL), nothing to do"
  exit 0
fi

IN=$(sed -n "${LINE}p" $LIST)
BASE=$(basename $IN)
NAME=${BASE%.fna.gz}
NAME=${NAME%.fna}
OUT=$OUTBASE/$NAME

echo "host   : $(hostname)"
echo "job    : $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID offset $OFFSET -> index $IDX"
echo "input  : $IN"
echo "name   : $NAME"
echo "start  : $(date)"

if [ -s "$OUT/${NAME}.json" ]; then
  echo "already complete, skipping"
  exit 0
fi

if [ ! -e "$IN" ]; then echo "INPUT MISSING: $IN"; exit 1; fi
if [ ! -x "$ASMASH" ]; then echo "ANTISMASH NOT EXECUTABLE: $ASMASH"; exit 1; fi

mkdir -p $OUTBASE $SCRATCH

FASTA=$IN
if [[ "$IN" == *.gz ]]; then
  FASTA=$SCRATCH/${NAME}.fna
  zcat $IN > $FASTA
  if [ ! -s $FASTA ]; then echo "DECOMPRESSION PRODUCED EMPTY FILE"; exit 1; fi
fi
echo "contigs: $(grep -c '^>' $FASTA)"

$ASMASH \
  --taxon bacteria \
  --genefinding-tool prodigal-m \
  --cb-knownclusters \
  --output-dir $OUT \
  --output-basename $NAME \
  --cpus 4 \
  $FASTA

RC=$?
echo "-----------------------------------------------------------"
echo "antismash exit: $RC"

if [[ "$IN" == *.gz ]]; then rm -f $FASTA; fi

if [ $RC -ne 0 ]; then echo "ANTISMASH FAILED for $NAME"; exit $RC; fi

JSON=$OUT/${NAME}.json
if [ ! -s "$JSON" ]; then echo "NO JSON OUTPUT for $NAME"; exit 1; fi

echo "end    : $(date)"
echo "regions: $(ls $OUT/*region*.gbk 2>/dev/null | wc -l)"
echo "ANTISMASH_FINISHED $NAME"
