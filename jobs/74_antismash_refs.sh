#!/bin/bash
# Reference genomes are searched with antiSMASH.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/46_antismash_refs.sh
# Output: work/bgc_refs/antismash/
# ANTISMASH_REFS_V2_20260805
#SBATCH --job-name=asmash_ref
#SBATCH --partition=stajichlab
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-999%60
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/asmash_ref_%A_%a.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/asmash_ref_%A_%a.err

# V1 (27245533) failed on all three test tasks with "Modules failing
# prerequisites": antismash could not find hmmpfam2, blastp or prodigal.
# Calling the antismash binary by absolute path is not enough because it
# shells out to helpers in the same env bin directory. The amphibian run
# (27242880/27242881) only worked because the submitting shell had
# `module load antismash/7.1.0` active and --export=ALL carried that PATH
# in, so that run was not reproducible as written. This sets PATH
# explicitly and checks every helper before starting.

ROOT=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
WORK=$ROOT/work/bgc_refs
LIST=$WORK/bgc_ref_input_list.txt
OUTBASE=$WORK/antismash
SCRATCH=$WORK/scratch
ASBIN=/opt/linux/rocky/8.x/x86_64/pkgs/antismash/7.1.0/bin
ASMASH=$ASBIN/antismash
export PATH=$ASBIN:$PATH

OFFSET=${OFFSET:-0}
IDX=$(( SLURM_ARRAY_TASK_ID + OFFSET ))
LINE=$(( IDX + 1 ))
TOTAL=$(wc -l < $LIST)
if [ $IDX -ge $TOTAL ]; then echo "index $IDX beyond $TOTAL"; exit 0; fi

IN=$(sed -n "${LINE}p" $LIST)
BASE=$(basename $IN)
NAME=${BASE%.fna.gz}; NAME=${NAME%.fna}
OUT=$OUTBASE/$NAME

echo "host   : $(hostname)"
echo "job    : $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID offset $OFFSET -> $IDX"
echo "input  : $IN"
echo "start  : $(date)"
for h in prodigal blastp hmmpfam2 hmmscan diamond; do
  P=$(command -v $h)
  echo "  helper $h: ${P:-NOT FOUND}"
  [ -z "$P" ] && { echo "HELPER MISSING ON PATH: $h"; exit 1; }
done

if [ -s "$OUT/${NAME}.json" ]; then echo "already complete, skipping"; exit 0; fi
if [ ! -e "$IN" ]; then echo "INPUT MISSING: $IN"; exit 1; fi

mkdir -p $OUTBASE $SCRATCH
FASTA=$IN
if [[ "$IN" == *.gz ]]; then
  FASTA=$SCRATCH/${NAME}.fna
  zcat $IN > $FASTA
  [ -s $FASTA ] || { echo "DECOMPRESSION EMPTY"; exit 1; }
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
echo "antismash exit: $RC"
if [[ "$IN" == *.gz ]]; then rm -f $FASTA; fi
if [ $RC -ne 0 ]; then echo "ANTISMASH FAILED for $NAME"; exit $RC; fi
if [ ! -s "$OUT/${NAME}.json" ]; then echo "NO JSON for $NAME"; exit 1; fi

echo "end    : $(date)"
echo "regions: $(ls $OUT/*region*.gbk 2>/dev/null | wc -l)"
echo "ANTISMASH_REF_FINISHED $NAME"
