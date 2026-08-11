#!/bin/bash
# Endotherm comparison genomes are searched with antiSMASH.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/52_antismash_endotherm.sh
# Output: work/bgc_endo/antismash/
# ANTISMASH_ENDOTHERM_V1_20260806
#SBATCH --job-name=asmash_endo
#SBATCH --partition=stajichlab
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-327%60
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/asmash_endo_%A_%a.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/asmash_endo_%A_%a.err

# THE ENDOTHERM GUT ARM R9 HAS NEVER HAD. Only two antiSMASH runs existed:
# 1,155 amphibian-derived genomes and 1,247 GTDB Ruminococcaceae references.
# GTDB is a TAXONOMIC comparator, not a host one (91.9% are themselves MAGs,
# isolation source mostly blank, see R7).
#
# 328 genomes = 280 EHI mammal/bird Ruminococcaceae + 48 Youngblut
# Ruminococcaceae, both selected on column 2 of their GTDB-Tk summaries.
#
# WHY EHI MAMMAL IS THE POINT: EHI mammals and EHI newts come from one
# consortium, one pipeline, one data release, so contiguity is matched BY
# CONSTRUCTION rather than by luck. That is the same logic that makes the
# 93.8% turnover control in R3 persuasive. It is the only host contrast in
# this dataset where the assembly-quality confound is controlled by design.
#
# WHY YOUNGBLUT IS SECONDARY: 79.2% of Youngblut Ruminococcaceae already sit
# inside a GTDB r220 species cluster, so the arm is not independent of the
# reference set, and 22 of 48 are chicken. Included because it is cheap.
# Do not build a claim on it alone.
#
# WHAT THIS TESTS: results/assembly_quality_arms.tsv shows complete-BGC
# recovery is a function of N50 that does NOT differ by catalog (at matched
# N50: 0.34/0.14/0.37 in the 20-40 kb bin, 0.67/0.65/0.68 in 40-80 kb). The
# question is whether it also does not differ by HOST once contiguity is
# matched by pipeline.
#
# PATH is set explicitly and every helper is checked. The amphibian run
# (27242880/27242881) only worked because the submitting shell had
# `module load antismash/7.1.0` active and --export=ALL carried it in.
# 328 genomes needs 328 tasks (0-327), one per genome, under the 1000 cap.

ROOT=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
WORK=$ROOT/work/bgc_endo
LIST=$WORK/bgc_endo_input_list.txt
OUTBASE=$WORK/antismash
SCRATCH=$WORK/scratch
ASBIN=/opt/linux/rocky/8.x/x86_64/pkgs/antismash/7.1.0/bin
ASMASH=$ASBIN/antismash
export PATH=$ASBIN:$PATH

OFFSET=${OFFSET:-0}
IDX=$(( SLURM_ARRAY_TASK_ID + OFFSET ))
LINE=$(( IDX + 1 ))

if [ ! -s "$LIST" ]; then echo "MISSING LIST: $LIST, run the staging script first"; exit 1; fi
TOTAL=$(wc -l < $LIST)
if [ $IDX -ge $TOTAL ]; then echo "index $IDX beyond $TOTAL"; exit 0; fi

IN=$(sed -n "${LINE}p" $LIST | cut -f1)
NAME=$(sed -n "${LINE}p" $LIST | cut -f2)
OUT=$OUTBASE/$NAME

echo "host   : $(hostname)"
echo "job    : $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID offset $OFFSET -> $IDX"
echo "input  : $IN"
echo "name   : $NAME"
echo "start  : $(date)"

for h in prodigal blastp hmmpfam2 hmmscan diamond; do
  P=$(command -v $h)
  echo "  helper $h: ${P:-NOT FOUND}"
  [ -z "$P" ] && { echo "HELPER MISSING ON PATH: $h"; exit 1; }
done

if [ -s "$OUT/${NAME}.json" ]; then echo "already complete, skipping"; exit 0; fi
if [ -d "$OUT" ]; then
  echo "output dir exists without a json, clearing it"
  rm -rf $OUT
fi
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
echo "ANTISMASH_ENDO_FINISHED $NAME"
