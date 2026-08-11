#!/bin/bash
# The joint tree is inferred with FastTree under LG. Gamma is not used, since it fails on this alignment.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/32_figure_tree.sh
# Output: work/rep_tree/figure_tree.nwk
# FIGURE_TREE_V1_20260804
#SBATCH --job-name=fig_tree
#SBATCH --partition=stajichlab
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/fig_tree.%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/fig_tree.%j.err

# Job 27143077 died with a PairLogLk numerical failure under -gamma on an
# alignment with large gap blocks. -gamma rescales branch lengths only and
# does not change topology, so it is dropped. Exit code captured directly:
# in 27143077 the reported exit came from the enclosing if-statement and a
# FastTree failure showed as exit 0.

ROOT=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
IN=$ROOT/work/rep_tree/figure_tree.faa
TMP=$ROOT/work/rep_tree/figure_tree.tmp.nwk
OUT=$ROOT/work/rep_tree/figure_tree.nwk
LOG=$ROOT/work/rep_tree/figure_tree_fasttree.log

module load fasttree
export OMP_NUM_THREADS=16

echo "host    : $(hostname)"
echo "job     : $SLURM_JOB_ID"
echo "started : $(date)"
echo "tips    : $(grep -c '^>' $IN)"
echo "-----------------------------------------------------------"

FastTreeMP -lg -log $LOG $IN > $TMP
RC=$?
echo "-----------------------------------------------------------"
echo "FastTree exit: $RC"
if [ $RC -ne 0 ]; then echo "FASTTREE FAILED"; exit $RC; fi
if [ ! -s $TMP ]; then echo "EMPTY TREE"; exit 1; fi
if ! tail -c 200 $TMP | grep -q ";"; then echo "NO TERMINATING SEMICOLON"; exit 1; fi
mv $TMP $OUT
echo "finished: $(date)"
echo "bytes   : $(stat -c%s $OUT)"
echo "FIGURE_TREE_FINISHED"
