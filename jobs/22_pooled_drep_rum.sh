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

# Reference completeness and contamination are taken from
# data/gtdb/bac120_metadata.tsv using checkm_completeness rather than CheckM2,
# matching how the MAG arms were scored, so one tool underlies the whole run.
# 1,224 of 1,247 references matched; 23 retain a stated fallback of
# completeness 95 and contamination 2.
# The -comp 70 floor drops 62 references before clustering. No MAG arm loses
# any genome. Cluster membership is therefore defined over the genomes that
# were compared, not over every genome staged.

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
NSTAGED=$(ls $GDIR/*.fa | wc -l)
echo "genomes staged : $NSTAGED"
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
echo "genomes clustered : $(( $(wc -l < $CDB) - 1 )) of $NSTAGED staged"
echo "secondary clusters: $(tail -n +2 $CDB | cut -d, -f2 | sort -u | wc -l)"
echo "Wdb rows          : $(( $(wc -l < $WDB) - 1 ))"
echo "POOLED_DREP_RUM_V3_FINISHED"
