#!/bin/bash
# Completeness and contamination are estimated for the focal genomes as the happi quality variable.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/35_checkm2_focal.sh
# Output: work/focal_genus_pangenome/checkm2_out/quality_report.tsv
# CHECKM2_FOCAL_V2_20260804
#SBATCH --job-name=ckm2_focal
#SBATCH --partition=stajichlab
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/ckm2_focal.%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/ckm2_focal.%j.err

# Uniform completeness across all 125 focal genomes with ONE tool. The
# amphibian arm currently carries CheckM1 values from the manifest and the
# reference arm carries GTDB metadata values, so any prevalence difference
# is confounded with how completeness was measured. happi needs one
# consistent quality covariate.
#
# The checkm2/1.1.0 module activates a conda env, and `module load` inside a
# batch shell produced "CondaError: Run 'conda init'" with exit 127 in job
# 27148950. The binary is therefore called by its resolved absolute path,
# which is the pattern that worked for dRep in job 27149754.

ROOT=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
WORK=$ROOT/work/focal_genus_pangenome
GDIR=$WORK/genomes
OUT=$WORK/checkm2_out

CHECKM2_BIN=/opt/linux/rocky/8.x/x86_64/pkgs/checkm2/1.1.0/bin/checkm2
export CHECKM2DB=/srv/projects/db/checkm2/CheckM2_database/uniref100.KO.1.dmnd

echo "host    : $(hostname)"
echo "job     : $SLURM_JOB_ID"
echo "started : $(date)"
echo "genomes : $(ls $GDIR/*.fna | wc -l)"
echo "binary  : $CHECKM2_BIN"
echo "db      : $CHECKM2DB"

if [ ! -x "$CHECKM2_BIN" ]; then echo "BINARY NOT EXECUTABLE"; exit 1; fi
if [ ! -e "$CHECKM2DB" ]; then
  echo "DB PATH NOT FOUND. Contents of the database directory:"
  ls -la /srv/projects/db/checkm2/CheckM2_database/
  exit 1
fi
$CHECKM2_BIN --version
echo "-----------------------------------------------------------"

mkdir -p $OUT

$CHECKM2_BIN predict \
  --input $GDIR \
  --extension .fna \
  --output-directory $OUT \
  --threads 24 \
  --force

RC=$?
echo "-----------------------------------------------------------"
echo "checkm2 exit: $RC"
if [ $RC -ne 0 ]; then echo "CHECKM2 FAILED"; exit $RC; fi

REP=$OUT/quality_report.tsv
if [ ! -s $REP ]; then echo "NO quality_report.tsv"; exit 1; fi

echo "finished : $(date)"
echo "genomes scored: $(( $(wc -l < $REP) - 1 )) of 125"
echo
echo "header:"
head -1 $REP
echo
echo "completeness and contamination by group (staged names are"
echo "GENUS__group__id, so field 2 of the underscore split is the group):"
awk -F'\t' 'NR>1 {
  split($1, a, "__"); g=a[2];
  cs[g]+=$2; xs[g]+=$3; n[g]++;
  if (n[g]==1 || $2<cmin[g]) cmin[g]=$2;
  if ($2>cmax[g]) cmax[g]=$2;
} END {
  for (g in n) printf "  %-12s n=%3d  completeness mean %.1f [%.1f, %.1f]  contamination mean %.2f\n", g, n[g], cs[g]/n[g], cmin[g], cmax[g], xs[g]/n[g];
}' $REP
echo
echo "If the amphibian and reference arms differ materially in completeness,"
echo "raw prevalence differences are confounded and happi must carry this as"
echo "the quality covariate."
echo "CHECKM2_FOCAL_FINISHED"
