#!/bin/bash
# Focal proteins are annotated with eggNOG-mapper for orthologous group assignment.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/38_eggnog_focal.sh
# Output: work/focal_genus_pangenome/eggnog/focal.emapper.annotations
# EGGNOG_FOCAL_V2_20260805
#SBATCH --job-name=eggnog
#SBATCH --partition=stajichlab
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/eggnog_focal.%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/eggnog_focal.%j.err

# Annotates CLUSTER REPRESENTATIVES, not all 334,627 proteins: every member
# of a cluster is the same gene family, so one annotation per cluster is
# enough and cuts the work roughly four-fold.
#
# id50 cov-mode 1 is the primary clustering. The merge test showed 81% of
# the reference-higher calls at id80 were split orthologs that rejoin at
# id50, and measured amphibian-to-reference identity is 62.7% (mean, p95
# 0.776), well below an 80% threshold, so id80 cannot cluster true
# orthologs across the two arms.
#
# The eggnog-mapper module activates a conda env; `module load` in a batch
# shell gave "CondaError: Run 'conda init'" and exit 127 in job 27148950.
# The binary is therefore called by absolute path and the database variable
# is set explicitly rather than inherited from the module.

ROOT=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
WORK=$ROOT/work/focal_genus_pangenome
FAA=$WORK/all_proteins.faa
CLU=$WORK/clu_covmode1_cluster.tsv
REPS=$WORK/cluster_reps_id50_cov1.faa
OUT=$WORK/eggnog

EMAPPER=/opt/linux/rocky/8.x/x86_64/pkgs/eggnog-mapper/2.1.9/env/bin/emapper.py
export EGGNOG_DATA_DIR=/srv/projects/db/eggNOG/LATEST
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/eggnog-mapper/2.1.9/env/bin:$PATH

echo "host  : $(hostname)"
echo "job   : $SLURM_JOB_ID"
echo "start : $(date)"
echo "db    : $EGGNOG_DATA_DIR"

if [ ! -f "$EMAPPER" ]; then echo "EMAPPER NOT FOUND: $EMAPPER"; exit 1; fi
if [ ! -d "$EGGNOG_DATA_DIR" ]; then echo "DB DIR NOT FOUND"; exit 1; fi
ls -la $EGGNOG_DATA_DIR | head -6
$EMAPPER --version
echo "-----------------------------------------------------------"

mkdir -p $OUT

cut -f1 $CLU | sort -u > $WORK/rep_ids.txt
echo "cluster representatives: $(wc -l < $WORK/rep_ids.txt)"

awk 'NR==FNR{keep[$1]; next}
     /^>/{n=substr($0,2); split(n,a," "); p=(a[1] in keep)}
     p' $WORK/rep_ids.txt $FAA > $REPS
echo "representative proteins written: $(grep -c '^>' $REPS)"
if [ ! -s $REPS ]; then echo "NO REPRESENTATIVES EXTRACTED"; exit 1; fi

$EMAPPER \
  -i $REPS \
  -o focal \
  --output_dir $OUT \
  -m diamond \
  --cpu 24 \
  --itype proteins \
  --tax_scope Bacteria \
  --go_evidence non-electronic \
  --target_orthologs all \
  --temp_dir $OUT \
  --override

RC=$?
echo "-----------------------------------------------------------"
echo "emapper exit: $RC"
if [ $RC -ne 0 ]; then echo "EGGNOG FAILED"; exit $RC; fi

ANN=$OUT/focal.emapper.annotations
if [ ! -s $ANN ]; then echo "NO ANNOTATIONS FILE"; exit 1; fi

echo "end   : $(date)"
echo "annotation rows: $(grep -vc '^#' $ANN)"
echo
echo "header (column order matters for parsing):"
grep '^#query' $ANN | head -1 | tr '\t' '\n' | cat -n
echo "EGGNOG_FOCAL_FINISHED"
