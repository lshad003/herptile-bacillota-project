#!/usr/bin/env bash
# DIAMOND blastp of amphibian cluster members against the reference Angelakisella proteomes.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/angelakisella_diamond.sh
# Output: work/angelakisella_diamond/hits.tsv
#SBATCH --job-name=angel_diamond
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/angelakisella_diamond.%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/angelakisella_diamond.%j.err

WORK=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/work/angelakisella_diamond
QUERY=$WORK/amphibian_queries.faa
REF=$WORK/reference_proteomes.faa
DB=$WORK/ref_db
HITS=$WORK/hits.tsv

module load diamond

DIAMOND=$(command -v diamond)
if [ -z "$DIAMOND" ]; then
    echo "STOP: diamond not found after module load" >&2
    exit 1
fi
echo "diamond binary: $DIAMOND"
$DIAMOND version

for f in "$QUERY" "$REF"; do
    if [ ! -s "$f" ]; then
        echo "STOP: missing or empty input: $f" >&2
        exit 1
    fi
done

if [ -s "$HITS" ]; then
    echo "STOP: hits file exists, refusing to overwrite: $HITS" >&2
    exit 1
fi

echo "query proteins: $(grep -c '^>' "$QUERY")"
echo "reference proteins: $(grep -c '^>' "$REF")"

$DIAMOND makedb --in "$REF" -d "$DB" --threads 8 || exit 1

$DIAMOND blastp --very-sensitive \
    -q "$QUERY" -d "$DB" -o "$HITS" \
    -e 1e-5 -k 5 --threads 8 \
    --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
    || exit 1

echo "hit rows: $(wc -l < "$HITS")"
echo "DONE"
# ANGELAKISELLA_DIAMOND_JOB_V1
