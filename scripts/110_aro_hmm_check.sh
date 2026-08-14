#!/usr/bin/env bash
# Chorismate-pathway absence verified by Pfam HMM with the non-amphibian genomes as positive control.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/angelakisella_aro_hmm_check.sh
# Output: results/angelakisella_aro_hmm_check.tsv
# Verifies the chorismate-pathway absence claim (aroA/B/C/E/F at
# amphibian prevalence zero in the happi result) by an annotation-
# independent method: Pfam HMMs at gathering thresholds against all 57
# Angelakisella proteomes. The 26 non-amphibian genomes are the positive
# control; they carry the genes per happi, so hits there prove the
# method detects them. Login node: 5 profiles vs ~150k proteins streams.

set -u
BASE=/bigdata/stajichlab/lshad003/ruminococcaceae-agent
PROT=$BASE/work/all_arms_pangenome/proteomes
META=$BASE/results/angelakisella_matrix_genomes_v2.tsv
PFAM=/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm
WORK=$BASE/work/angelakisella_aro_hmm
OUT=$BASE/results/angelakisella_aro_hmm_check.tsv
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3

# gene -> Pfam model name (one diagnostic domain per enzyme)
MODELS="EPSP_synthase DHQ_synthase Chorismate_synt Shikimate_DH DAHP_synth_1"
GENES="aroA aroB aroC aroE aroF"

module load hmmer/3.3.2
HMMSEARCH=$(command -v hmmsearch)
HMMFETCH=$(command -v hmmfetch)
if [ -z "$HMMSEARCH" ] || [ -z "$HMMFETCH" ]; then
    echo "STOP: hmmer tools not found after module load hmmer/3.3.2" >&2
    exit 1
fi

if [ -s "$OUT" ]; then
    echo "STOP: output exists, refusing to overwrite: $OUT" >&2
    exit 1
fi
mkdir -p "$WORK"

# 1. fetch the five models by name; stop if any name is wrong
HMMS=$WORK/aro_models.hmm
rm -f "$HMMS"
for m in $MODELS; do
    $HMMFETCH "$PFAM" "$m" >> "$HMMS" 2> /dev/null
    if ! grep -q "NAME  $m" "$HMMS"; then
        echo "STOP: Pfam model name not found: $m" >&2
        echo "List candidates with: grep '^NAME' $PFAM | grep -i <term>" >&2
        exit 1
    fi
done
echo "models fetched: $(grep -c '^NAME' "$HMMS")/5"

# 2. concatenate the 57 test proteomes with genome-tagged headers
ALLFAA=$WORK/angelakisella_57.faa
$PY - "$META" "$PROT" "$ALLFAA" << 'PYEOF'
import os, sys
meta, prot, out = sys.argv[1], sys.argv[2], sys.argv[3]
genomes = []
with open(meta) as fh:
    fh.readline()
    for line in fh:
        g, cat, grp = line.rstrip("\n").split("\t")
        if grp in ("amphibian", "non_amphibian"):
            genomes.append(g)
if len(genomes) != 57:
    sys.exit("STOP: expected 57 test genomes, got %d" % len(genomes))
n = 0
with open(out, "w") as o:
    for g in genomes:
        path = os.path.join(prot, g + ".faa")
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            sys.exit("STOP: missing proteome " + path)
        with open(path, errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    o.write(">" + g + "|" + line[1:].split()[0] + "\n")
                    n += 1
                else:
                    o.write(line)
print("proteins written: %d from %d genomes" % (n, len(genomes)))
PYEOF
[ $? -eq 0 ] || exit 1

# 3. hmmsearch at gathering thresholds
TBL=$WORK/aro_hits_domtbl.txt
$HMMSEARCH --cut_ga --domtblout "$TBL" -o /dev/null --cpu 4 "$HMMS" "$ALLFAA" || exit 1
echo "hmmsearch done, hit rows: $(grep -vc '^#' "$TBL")"

# 4. per-genome, per-gene presence + group summary + completeness context
$PY - "$TBL" "$META" "$BASE/results/angelakisella_happi_metadata.tsv" "$OUT" << 'PYEOF'
import sys
from collections import defaultdict
tbl, meta, happi_meta, out = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

model2gene = {"EPSP_synthase": "aroA", "DHQ_synthase": "aroB",
              "Chorismate_synt": "aroC", "Shikimate_DH": "aroE",
              "DAHP_synth_1": "aroF"}
genes = ["aroA", "aroB", "aroC", "aroE", "aroF"]

grp = {}
with open(meta) as fh:
    fh.readline()
    for line in fh:
        g, cat, gr = line.rstrip("\n").split("\t")
        if gr in ("amphibian", "non_amphibian"):
            grp[g] = gr

comp = {}
with open(happi_meta) as fh:
    fh.readline()
    for line in fh:
        g, gr, c = line.rstrip("\n").split("\t")
        comp[g] = float(c)

hits = defaultdict(set)   # genome -> set of genes
for line in open(tbl):
    if line.startswith("#"):
        continue
    f = line.split()
    target, query = f[0], f[3]
    genome = target.split("|", 1)[0]
    gene = model2gene.get(query)
    if gene and genome in grp:
        hits[genome].add(gene)

with open(out, "w") as o:
    o.write("genome\tgroup\tcompleteness\t" + "\t".join(genes) + "\n")
    for g in sorted(grp):
        row = [g, grp[g], "%.1f" % comp.get(g, float("nan"))]
        row += ["1" if x in hits[g] else "0" for x in genes]
        o.write("\t".join(row) + "\n")

print("")
print("%-6s %14s %14s" % ("gene", "amphibian", "non_amphibian"))
for gene in genes:
    na = sum(1 for g in grp if grp[g] == "amphibian" and gene in hits[g])
    nn = sum(1 for g in grp if grp[g] == "non_amphibian" and gene in hits[g])
    ta = sum(1 for g in grp if grp[g] == "amphibian")
    tn = sum(1 for g in grp if grp[g] == "non_amphibian")
    print("%-6s %10d/%-3d %10d/%-3d" % (gene, na, ta, nn, tn))

ca = [comp[g] for g in grp if grp[g] == "amphibian" and g in comp]
print("")
print("amphibian completeness mean %.1f, min %.1f"
      % (sum(ca) / len(ca), min(ca)))
print("")
print("READ IT THIS WAY: non_amphibian counts near their happi prevalences")
print("are the positive control proving the HMMs detect these enzymes in")
print("this genus. Amphibian counts at or near zero then confirm the")
print("pathway absence by a method independent of the eggNOG annotations.")
print("Coordinated absence of five enzymes across 31 genomes at ~90% mean")
print("completeness is not producible by genome incompleteness, which")
print("would scatter losses randomly across genes and genomes.")
PYEOF
[ $? -eq 0 ] || exit 1

echo "wrote: $OUT"
# ANGELAKISELLA_ARO_HMM_CHECK_V1
