#!/bin/bash
# The five staged arms are dereplicated together at 95 percent ANI with dRep 3.5.0.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/30_pooled_drep_rum.sh
# Output: work/pooled_drep_rum/drep_out/data_tables/Cdb.csv and Wdb.csv
# POOLED_DREP_RUM_V3_20260804
#SBATCH --job-name=drep_rum
#SBATCH --partition=stajichlab
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/drep_rum.%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/drep_rum.%j.err

# V1 (27148950) failed: `module load drep` gives CondaError in a batch shell,
#   exit 127. Environment now copied verbatim from scripts/drep_rerun_2229.sh.
# V2 (27149185) pended indefinitely: epyc fully allocated. Cancelled.
# V3: stajichlab partition (r11 idle), plus 1,247 GTDB r220 Ruminococcaceae
#   references staged in, so 1,795 genomes total.
#   Reference completeness/contamination taken from data/gtdb/bac120_metadata.tsv
#   (checkm_completeness, NOT checkm2, matching how the MAGs were scored).
#   1,224 of 1,247 matched; 23 retain a stated fallback of comp 95 / con 2.
#   62 references fall below the -comp 70 floor and will be dropped, so
#   expect 1,733 genomes clustered from 1,795 staged. No MAG arm loses any.

BASE=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
DREPENV=/bigdata/stajichlab/lshad003/condaenvs/drep
DREPBIN=$DREPENV/bin
WORK=$BASE/work/pooled_drep_rum
GDIR=$WORK/genomes
INFO=$WORK/genome_info.csv
OUT=$WORK/drep_out
export PATH=$DREPBIN:$PATH

echo "host    : $(hostname)"
echo "job     : $SLURM_JOB_ID"
echo "started : $(date)"
echo
echo "=== versions ==="
ls -1 $DREPENV/conda-meta/ | grep -i -E "^drep|^fastani|^mash"
echo
$DREPBIN/dRep check_dependencies
echo
echo "genomes staged : $(ls $GDIR/*.fa | wc -l)"
echo "genomeInfo rows: $(( $(wc -l < $INFO) - 1 ))"
echo "-----------------------------------------------------------"

# Parameters identical to scripts/drep_rerun_2229.sh, the run that produced
# the 1,171 herptile SGBs, so the pooled clusters are directly comparable.

$DREPBIN/dRep dereplicate $OUT \
  -g $GDIR/*.fa \
  --genomeInfo $INFO \
  -comp 70 -con 10 -l 50000 \
  -pa 0.90 -sa 0.95 -nc 0.30 --coverage_method larger \
  --S_algorithm fastANI --clusterAlg average \
  -p 24

RC=$?
echo "-----------------------------------------------------------"
echo "dRep exit code: $RC"
if [ $RC -ne 0 ]; then echo "DREP FAILED"; exit $RC; fi

CDB=$OUT/data_tables/Cdb.csv
WDB=$OUT/data_tables/Wdb.csv
if [ ! -s $CDB ]; then echo "NO Cdb.csv PRODUCED"; exit 1; fi

echo "finished : $(date)"
echo
ls -1 $OUT/data_tables
echo
echo "genomes clustered : $(( $(wc -l < $CDB) - 1 )) of 1795 staged"
echo "secondary clusters: $(tail -n +2 $CDB | cut -d, -f2 | sort -u | wc -l)"
echo "Wdb rows          : $(( $(wc -l < $WDB) - 1 ))"
echo "POOLED_DREP_RUM_V3_FINISHED"
