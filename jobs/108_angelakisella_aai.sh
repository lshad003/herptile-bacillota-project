#!/bin/bash
# Amphibian-to-reference amino acid identity measured over 100 genome pairs.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/angelakisella_aai.sh
# Output: results/angelakisella_aai_pairs.tsv
#SBATCH --job-name=angel_aai
#SBATCH --partition=short
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/angelakisella_aai.%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/angelakisella_aai.%j.err

# Genome-wide amphibian-vs-reference AAI for Angelakisella, so the
# over-splitting result is reported at a measured divergence rather than
# a value borrowed from the focal genera. 10 amphibian x 10 reference
# genomes, two-way DIAMOND best hits, mean identity of reciprocal pairs.
# sbatch per the DIAMOND rule; short partition, well under the cap.

set -u
BASE=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
PROT=$BASE/work/all_arms_pangenome/proteomes
WORK=$BASE/work/angelakisella_aai
OUT=$BASE/results/angelakisella_aai_pairs.tsv
DIAMOND=/opt/linux/rocky/8.x/x86_64/pkgs/diamond/2.1.24/bin/diamond
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3

if [ ! -x "$DIAMOND" ]; then echo "STOP: diamond not found" >&2; exit 1; fi
if [ -s "$OUT" ]; then echo "STOP: output exists: $OUT" >&2; exit 1; fi
mkdir -p "$WORK"

# 10 amphibian (5 herptile + 5 newt, alphabetical for reproducibility)
AMPH="herptile__STP248.12601_R.bin.113 herptile__STP248.12601_R.bin.126 herptile__STP572.12604_R.bin.29 herptile__UHM1073.41096_R.bin.113 herptile__UHM1088.23067_R.bin.62 ehi_amphibian__EHM034541 ehi_amphibian__EHM034833 ehi_amphibian__EHM035004 ehi_amphibian__EHM047839 ehi_amphibian__EHM047946"
REF="gtdb_ref__GCA_003453215.1 gtdb_ref__GCA_004554485.1 gtdb_ref__GCA_004557855.1 gtdb_ref__GCA_013316495.1 gtdb_ref__GCA_019423445.1 gtdb_ref__GCA_020024515.1 gtdb_ref__GCA_020024555.1 gtdb_ref__GCA_020025195.1 gtdb_ref__GCA_021769525.1 gtdb_ref__GCA_022794295.1"

for g in $AMPH $REF; do
    if [ ! -s "$PROT/$g.faa" ]; then
        echo "STOP: missing proteome $PROT/$g.faa" >&2; exit 1
    fi
done

for g in $REF; do
    if [ ! -s "$WORK/$g.dmnd" ]; then
        $DIAMOND makedb --in "$PROT/$g.faa" -d "$WORK/$g" --threads 8 --quiet || exit 1
    fi
done
for g in $AMPH; do
    if [ ! -s "$WORK/$g.dmnd" ]; then
        $DIAMOND makedb --in "$PROT/$g.faa" -d "$WORK/$g" --threads 8 --quiet || exit 1
    fi
done

echo "pair AAI (two-way best hits, qcov>=50, id>=20, e<=1e-5):"
echo -e "amphibian\treference\tn_rbh\taai" > "$OUT"
for a in $AMPH; do
  for r in $REF; do
    F1=$WORK/f_${a}__VS__${r}.tsv
    F2=$WORK/r_${a}__VS__${r}.tsv
    $DIAMOND blastp --very-sensitive --quiet -q "$PROT/$a.faa" -d "$WORK/$r" \
      -o "$F1" -e 1e-5 -k 1 --threads 8 \
      --outfmt 6 qseqid sseqid pident length qlen || exit 1
    $DIAMOND blastp --very-sensitive --quiet -q "$PROT/$r.faa" -d "$WORK/$a" \
      -o "$F2" -e 1e-5 -k 1 --threads 8 \
      --outfmt 6 qseqid sseqid pident length qlen || exit 1
    $PY - "$F1" "$F2" "$a" "$r" "$OUT" << 'PYEOF'
import sys
f1, f2, a, r, out = sys.argv[1:6]
def best(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            q, s, pid, alen, qlen = line.rstrip("\n").split("\t")
            if 100.0 * float(alen) / float(qlen) < 50.0:
                continue
            if float(pid) < 20.0:
                continue
            if q not in d:
                d[q] = (s, float(pid))
    return d
fwd = best(f1)
rev = best(f2)
idents = []
for q, (s, pid) in fwd.items():
    z = rev.get(s)
    if z is not None and z[0] == q:
        idents.append((pid + z[1]) / 2.0)
aai = sum(idents) / len(idents) if idents else float("nan")
with open(out, "a") as o:
    o.write("%s\t%s\t%d\t%.2f\n" % (a, r, len(idents), aai))
print("%s vs %s: n_rbh=%d aai=%.2f" % (a, r, len(idents), aai))
PYEOF
    rm -f "$F1" "$F2"
  done
done

echo "---"
$PY - "$OUT" << 'PYEOF'
import sys
vals = []
with open(sys.argv[1]) as fh:
    fh.readline()
    for line in fh:
        vals.append(float(line.rstrip("\n").split("\t")[3]))
vals.sort()
n = len(vals)
print("pairs: %d" % n)
print("AAI mean %.2f, median %.2f, min %.2f, max %.2f"
      % (sum(vals) / n, vals[n // 2], vals[0], vals[-1]))
PYEOF
echo "DONE"
# ANGELAKISELLA_AAI_JOB_V1
