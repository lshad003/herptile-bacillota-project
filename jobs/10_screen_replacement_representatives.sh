#!/bin/bash -l
# Classification and screening of two replacement SGB representatives.
#
# Source: ruminococcaceae-agent/jobs/67_repick_two_reps.sh
# Writes: results/gtdbtk_repick_two_r220/, results/gunc_repick_two/,
#         results/rrna_trna_repick_two/
#
# SGBs 351_1 and 382_2 contain genomes from wild-caught animals but their
# dRep representatives were drawn from the laboratory-reared arm. When that
# arm is excluded, both require a new representative from the remaining
# members. Neither replacement had previously been a representative, so
# neither had been classified or screened.
#
# CheckM is not rerun. Completeness and contamination for all MAGs come from
# a single CheckM1 run recorded in the MAG manifest; running a different tool
# on two genomes would make their quality values incomparable with the rest.
# REPICK_TWO_REPS_V1_20260808
#SBATCH --job-name=repick2
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/repick2_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/repick2_%j.err
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=06:00:00

# SBATCH RATIONALE: GTDB-Tk classify_wf loads the full r220 reference data
# into memory regardless of genome count, so the 500 G and highmem from job
# 27138736 are kept. That job did 718 genomes in 1h21m with these settings;
# 6 h is margin for two. GUNC's DIAMOND against the 13.5 GB proGenomes db is
# the other memory consumer and runs in the same allocation.
#
# PURPOSE: SGBs 351_1 and 382_2 are wild-containing but their dRep
# representatives are VIVARIUM (wood frog) MAGs. If the wood frog arm is
# excluded from the paper, both need a new representative from the remaining
# members. Neither replacement has ever been a representative, so neither has
# been through GTDB-Tk, GUNC, barrnap or tRNAscan.
#
# CheckM is deliberately NOT re-run. The manifest already carries CheckM1
# completeness and contamination for all 2,229 MAGs including these two, from
# the same run as every other genome. Running CheckM2 on two genomes would
# mix tools within one comparison.
#
# Log paths ABSOLUTE (27102139 died with no log on a relative one).
# No `module purge` (breaks the module system in batch shells).
# No `set -euo pipefail`.
# DO NOT `module load diamond`: GUNC 1.0.6 pins DIAMOND to 2.0.4 and ships it.

ROOT=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
WORK=$ROOT/work/repick_two
OUT_GTDB=$ROOT/results/gtdbtk_repick_two_r220
OUT_GUNC=$ROOT/results/gunc_repick_two
OUT_FEAT=$ROOT/results/rrna_trna_repick_two
BATCH=$WORK/batchfile.tsv
GUNCIN=$WORK/gunc_input
TMPG=$WORK/tmp_gtdbtk
TMPU=$WORK/tmp_gunc

A=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results/UHM245.23042/bins/UHM245.23042_R.bin.35.fa
B=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results/UHM975.23062/bins/UHM975.23062_R.bin.174.fa

PRODIGAL=/opt/linux/rocky/8.x/x86_64/pkgs/antismash/7.1.0/bin/prodigal
BARRNAP=/opt/linux/rocky/8.x/x86_64/pkgs/barrnap/0.9/bin/barrnap
TRNASCAN=/opt/linux/rocky/8.x/x86_64/pkgs/trnascan-se/2.0.12/bin/tRNAscan-SE
DB=/srv/projects/db/GUNC/gunc_db_progenomes2.1.dmnd

echo "[START] $(date) host=$(hostname) job=$SLURM_JOB_ID"

for f in "$A" "$B"; do
  if [ ! -s "$f" ]; then echo "[FATAL] missing genome: $f" >&2; exit 1; fi
done
for b in "$PRODIGAL" "$BARRNAP" "$TRNASCAN"; do
  if [ ! -x "$b" ]; then echo "[FATAL] not executable: $b" >&2; exit 1; fi
done
if [ ! -s "$DB" ]; then echo "[FATAL] GUNC database missing: $DB" >&2; exit 1; fi

PV=$($PRODIGAL -v 2>&1 | grep -oE 'V[0-9]+\.[0-9]+\.[0-9]+' | head -1)
echo "[INFO] prodigal $PV"
if [ "$PV" != "V2.6.3" ]; then echo "[FATAL] expected Prodigal V2.6.3, found $PV" >&2; exit 1; fi

rm -rf $WORK
mkdir -p $WORK $GUNCIN $TMPG $TMPU $OUT_GTDB $OUT_GUNC $OUT_FEAT

printf "%s\t%s\n" "$A" "UHM245.23042_R.bin.35"   >  $BATCH
printf "%s\t%s\n" "$B" "UHM975.23062_R.bin.174"  >> $BATCH
echo "[INFO] batchfile rows: $(wc -l < $BATCH)"
if [ "$(wc -l < $BATCH)" -ne 2 ]; then echo "[FATAL] batchfile not 2 rows" >&2; exit 1; fi

# ---------------------------------------------------------- GTDB-Tk
module load gtdbtk/2.4.1
export GTDBTK_DATA_PATH=/srv/projects/db/gtdbtk/220
echo "[INFO] gtdbtk: $(which gtdbtk)"
gtdbtk --version

gtdbtk classify_wf \
  --batchfile $BATCH \
  --out_dir $OUT_GTDB \
  --cpus 16 \
  --pplacer_cpus 8 \
  --skip_ani_screen \
  --tmpdir $TMPG
echo "[INFO] gtdbtk rc=$?"

SUMG=$OUT_GTDB/gtdbtk.bac120.summary.tsv
if [ -s "$SUMG" ]; then
  echo "[INFO] gtdbtk rows: $(( $(wc -l < $SUMG) - 1 )) of 2"
else
  echo "[WARN] no gtdbtk summary at $SUMG"
fi

# ---------------------------------------------------------- Prodigal for GUNC
$PRODIGAL -i $A -a $GUNCIN/UHM245.23042_R.bin.35.faa   -p meta -q > /dev/null 2>&1
$PRODIGAL -i $B -a $GUNCIN/UHM975.23062_R.bin.174.faa  -p meta -q > /dev/null 2>&1
NP=$(ls -1 $GUNCIN/*.faa 2>/dev/null | wc -l)
echo "[INFO] proteomes called: $NP"
if [ "$NP" -ne 2 ]; then echo "[FATAL] expected 2 proteomes, got $NP" >&2; exit 1; fi
for f in $GUNCIN/*.faa; do
  echo "[INFO] $(basename $f): $(grep -c '^>' $f) proteins"
done

# ---------------------------------------------------------- GUNC
source /etc/profile.d/modules.sh || true
module load gunc
if ! command -v gunc > /dev/null; then echo "[FATAL] gunc not on PATH" >&2; exit 1; fi
DV=$(diamond --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
echo "[INFO] gunc: $(command -v gunc), diamond $DV"
if [ "$DV" != "2.0.4" ]; then echo "[FATAL] GUNC needs diamond 2.0.4, found $DV" >&2; exit 1; fi

export TMPDIR=$TMPU
gunc run \
  --input_dir $GUNCIN \
  --gene_calls \
  --file_suffix .faa \
  --db_file $DB \
  --out_dir $OUT_GUNC \
  --threads ${SLURM_CPUS_PER_TASK} \
  --temp_dir $TMPU \
  --detailed_output
echo "[INFO] gunc rc=$?"

SUMU=$(ls -1 $OUT_GUNC/GUNC.*.maxCSS_level.tsv 2>/dev/null | head -1)
if [ -s "$SUMU" ]; then
  echo "[INFO] gunc summary: $SUMU"
  cat $SUMU
else
  echo "[WARN] no maxCSS_level.tsv in $OUT_GUNC"
fi

# ---------------------------------------------------------- barrnap + tRNAscan
mkdir -p $OUT_FEAT/barrnap $OUT_FEAT/trnascan
for f in "$A" "$B"; do
  ID=$(basename $f .fa)
  $BARRNAP --kingdom bac --quiet "$f" > $OUT_FEAT/barrnap/${ID}.gff 2>/dev/null
  echo "[INFO] barrnap $ID: $(grep -vc '^#' $OUT_FEAT/barrnap/${ID}.gff) feature lines"
  $TRNASCAN -B -o $OUT_FEAT/trnascan/${ID}.txt -q "$f" > /dev/null 2>&1
  if [ -s "$OUT_FEAT/trnascan/${ID}.txt" ]; then
    echo "[INFO] tRNAscan $ID: $(( $(wc -l < $OUT_FEAT/trnascan/${ID}.txt) - 3 )) tRNAs"
  else
    echo "[WARN] tRNAscan produced nothing for $ID"
  fi
done

echo "-----------------------------------------------------------"
echo "[DONE] $(date)"
echo "  gtdbtk : $OUT_GTDB"
echo "  gunc   : $OUT_GUNC"
echo "  rrna   : $OUT_FEAT"
echo "REPICK_TWO_REPS_FINISHED"
