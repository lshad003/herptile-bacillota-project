#!/bin/bash -l
# Dereplication of 2,229 Bacillota_A MAGs into species-level genome bins.
#
# Source: ruminococcaceae-agent/scripts/drep_rerun_2229.sh
# Run:    2026-07-23
# Output: results/drep_herptile_95ani_2229/
#
# Parameters were confirmed against dRep's own record of the executed run,
# results/drep_herptile_95ani_2229/log/cluster_arguments.json. An earlier
# script, jobs/14_drep_herptile.sh, writes to a different output directory
# and does not correspond to the run reported in the manuscript.
#
# Representatives are selected by dRep using completeness weight 1,
# contamination weight 5, N50 weight 0.5 and centrality weight 1.
#SBATCH --job-name=drep2229
#SBATCH --partition=epyc
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/drep_rerun_2229.%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/drep_rerun_2229.%j.err

BASE=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
DREPENV=/bigdata/stajichlab/lshad003/condaenvs/drep
DREPBIN=$DREPENV/bin
OUT=$BASE/results/drep_herptile_95ani_2229

export PATH=$DREPBIN:$PATH

echo "=== versions ==="
ls -1 $DREPENV/conda-meta/ | grep -i -E "^drep|^fastani|^mash"
echo
$DREPBIN/dRep check_dependencies
echo
echo "genomes in list: $(wc -w < $BASE/data/drep_rerun_2229_genome_list.txt)"
echo

$DREPBIN/dRep dereplicate $OUT \
  -g $(cat $BASE/data/drep_rerun_2229_genome_list.txt) \
  --genomeInfo $BASE/data/checkm_for_drep_2229.csv \
  -comp 70 -con 10 -l 50000 \
  -pa 0.90 -sa 0.95 -nc 0.30 --coverage_method larger \
  --S_algorithm fastANI --clusterAlg average \
  -p 24

echo "EXIT: $?"
echo
ls -1 $OUT/data_tables
echo
echo "SGB count (Wdb rows minus header):"
wc -l $OUT/data_tables/Wdb.csv
