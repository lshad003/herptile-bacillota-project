#!/usr/bin/env bash
# The absence of the family from reference genomes is closed by tblastn against the reference assemblies.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/jobs/angelakisella_nohit_tblastn.sh
# Output: work/angelakisella_diamond/nohit_probe_tblastn.tsv
#SBATCH --job-name=angel_tblastn
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/angelakisella_nohit_tblastn.%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ruminococcaceae-agent/logs/angelakisella_nohit_tblastn.%j.err

# Closes the last loophole on the absence claim for the 532-member family:
# tblastn of the 597 aa amphibian probe against the nucleotide contigs of
# the 25 reference Angelakisella genomes. Detects gene copies Prodigal
# missed (contig edges, frameshifts, unpredicted ORFs).
# V2: gtdb_ref genome files are gzipped (.fa.gz), verified by ls.

BASE=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
WORK=$BASE/work/angelakisella_diamond
GENOMES=$BASE/work/all_arms_pangenome/genomes
META=$BASE/work/rep_tree/figure_tree_metadata_genus.tsv
PROBE=$WORK/nohit_probe.faa
REFNUC=$WORK/ref_genomes_25.fna
DB=$WORK/ref_genomes_25_db
OUT=$WORK/nohit_probe_tblastn.tsv
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3

module load ncbi-blast

TBLASTN=$(command -v tblastn)
MAKEBLASTDB=$(command -v makeblastdb)
if [ -z "$TBLASTN" ] || [ -z "$MAKEBLASTDB" ]; then
    echo "STOP: tblastn/makeblastdb not found after module load ncbi-blast." >&2
    exit 1
fi
echo "tblastn: $TBLASTN"
$TBLASTN -version

if [ ! -s "$PROBE" ]; then
    echo "STOP: probe FASTA missing: $PROBE" >&2
    exit 1
fi

if [ -s "$OUT" ]; then
    echo "STOP: output exists, refusing to overwrite: $OUT" >&2
    exit 1
fi

$PY - "$META" "$GENOMES" "$REFNUC" << 'PYEOF'
import csv
import gzip
import os
import sys

meta, gdir, out = sys.argv[1], sys.argv[2], sys.argv[3]

KEEP_ARMS = {"herptile", "ehi_amphibian", "gtdb_ref"}


def core_genome_id(x):
    x = x.strip()
    if "__" in x:
        prefix, rest = x.split("__", 1)
        if prefix in KEEP_ARMS:
            x = rest
    elif "|" in x:
        prefix, rest = x.split("|", 1)
        if prefix in KEEP_ARMS:
            x = rest
    if "|" in x:
        x = x.split("|", 1)[0]
    if x.startswith(("RS_", "GB_")):
        x = x[3:]
    return x


refs = set()
with open(meta) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["genus"].strip() != "Angelakisella":
            continue
        if r["arm"].strip() != "gtdb_ref":
            continue
        refs.add(core_genome_id(r["genome"]))

print("reference Angelakisella genomes in metadata:", len(refs))
if len(refs) != 25:
    sys.exit("STOP: expected 25, got %d" % len(refs))

stems = {}
for entry in os.scandir(gdir):
    name = entry.name
    if name.endswith(".fa.gz"):
        stem = name[:-6]
        gzipped = True
    elif name.endswith(".fa"):
        stem = name[:-3]
        gzipped = False
    else:
        continue
    if not stem.startswith("gtdb_ref__"):
        continue
    gid = core_genome_id(stem.split("__", 1)[1])
    if gid in refs:
        stems[stem] = (entry.path, gzipped)

print("matching genome FASTA files found:", len(stems))
if len(stems) != 25:
    found = {core_genome_id(s.split("__", 1)[1]) for s in stems}
    for miss in sorted(refs - found):
        print("  missing:", miss)
    sys.exit("STOP: expected 25 genome files, got %d" % len(stems))

n = 0
with open(out, "w") as o:
    for stem in sorted(stems):
        path, gzipped = stems[stem]
        opener = gzip.open if gzipped else open
        with opener(path, "rt", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    o.write(">" + stem + "|" + line[1:].split()[0] + "\n")
                    n += 1
                else:
                    o.write(line)

print("contigs written:", n)
PYEOF
[ $? -eq 0 ] || exit 1

$MAKEBLASTDB -in "$REFNUC" -dbtype nucl -out "$DB" || exit 1

$TBLASTN -query "$PROBE" -db "$DB" -out "$OUT" \
    -evalue 10 -num_threads 4 \
    -outfmt "6 qseqid sseqid pident length qlen sstart send evalue bitscore" \
    || exit 1

echo "---"
echo "tblastn rows: $(wc -l < "$OUT")"
cat "$OUT"
echo "---"
echo "Ceiling is e=10, so weak junk WILL appear. What matters: any hit with"
echo "high identity over a long stretch means a reference genome carries the"
echo "gene unannotated and the absence claim dies. Only short low-identity"
echo "scraps means the absence holds at the DNA level in all 25 genomes."
echo "DONE"
# ANGELAKISELLA_NOHIT_TBLASTN_JOB_V2
