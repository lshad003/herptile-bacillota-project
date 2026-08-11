# Genera are tested for amphibian and reference interleaving by Fitch parsimony with 499 permutations.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/check_genus_interleaving.py
# Output: results/genus_interleaving.tsv
# CHECK_GENUS_INTERLEAVING_V3_20260804
import os, sys, random
from collections import Counter, defaultdict

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
TREE = os.path.join(ROOT, "work/rep_tree/figure_tree.nwk")
META = os.path.join(ROOT, "work/rep_tree/figure_tree_metadata.tsv")
GTDB_TAX = "/srv/projects/db/gtdbtk/220/taxonomy/gtdb_taxonomy.tsv"
OUT = os.path.join(ROOT, "results/genus_interleaving.tsv")

N_PERM = 499
random.seed(20260804)
AMPH_ARMS = {"herptile", "ehi_amphibian"}
MIN_PER_SIDE = 3

import dendropy
print("dendropy %s" % dendropy.__version__)


def read_tsv(path):
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                f = f + [""] * (len(header) - len(f))
            rows.append(dict(zip(header, f)))
    return rows


def parse_tax(s):
    out = {}
    for part in s.strip().split(";"):
        part = part.strip()
        if len(part) > 3 and part[1:3] == "__":
            out[part[0]] = part[3:]
    return out


def fitch(t2, state):
    """Fitch parsimony score. V1/V2 counted edges where the child and parent
    state sets were disjoint; because Fitch assigns the union when children
    disagree, and a union always intersects each child, that test never
    fired and the score was always 0. The score is the number of union
    events during the up-pass. Verified against four sanity cases."""
    up = {}
    score = 0
    for nd in t2.postorder_node_iter():
        if nd.is_leaf():
            lab = nd.taxon.label if nd.taxon is not None else None
            up[nd] = {state[lab]} if lab in state else None
            continue
        kids = [up[c] for c in nd.child_nodes() if up.get(c) is not None]
        if not kids:
            up[nd] = None
        elif len(kids) == 1:
            up[nd] = set(kids[0])
        else:
            inter = set.intersection(*kids)
            if inter:
                up[nd] = inter
            else:
                up[nd] = set.union(*kids)
                score += 1
    return score


meta = {r["tip"]: r for r in read_tsv(META)}

# The reference tips are labelled gtdb_ref|gtdb_ref|ACCESSION: the FASTA
# headers written by extract_ref_msa_rum.py already carried a gtdb_ref
# prefix and build_figure_tree_msa.py prefixed them a second time. Take
# the LAST pipe field as the accession.
ref_genus = {}
with open(GTDB_TAX) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        t = parse_tax(f[1])
        if t.get("f", "") == "Ruminococcaceae":
            acc = f[0].strip().replace("GB_", "").replace("RS_", "")
            ref_genus[acc] = t.get("g", "")

n_fixed = 0
for tip, m in meta.items():
    if m["arm"] != "gtdb_ref":
        continue
    acc = tip.split("|")[-1]
    g = ref_genus.get(acc, "")
    if g:
        m["genus"] = g
        n_fixed += 1
print("reference tips given a genus label: %d" % n_fixed)
if n_fixed == 0:
    sys.exit("reference genus join still failing, stopping")

tree = dendropy.Tree.get(path=TREE, schema="newick",
                         preserve_underscores=True, rooting="force-unrooted")
tips = [t.label for t in tree.taxon_namespace]
print("tree tips: %d" % len(tips))


def group_of(m):
    if m["arm"] in AMPH_ARMS:
        return "amphibian"
    if m["arm"] == "gtdb_ref":
        return "reference"
    if m["arm"] == "outgroup":
        return "outgroup"
    return "endotherm"


by_genus = defaultdict(lambda: defaultdict(list))
for t in tips:
    m = meta.get(t)
    if not m or not m["genus"] or m["arm"] == "outgroup":
        continue
    by_genus[m["genus"]][group_of(m)].append(t)

testable = []
for g, d in by_genus.items():
    a = len(d.get("amphibian", []))
    o = len(d.get("endotherm", [])) + len(d.get("reference", []))
    if a >= MIN_PER_SIDE and o >= MIN_PER_SIDE:
        testable.append((g, a, o))
testable.sort(key=lambda x: -x[1])

print()
print("=" * 74)
print("GENERA WITH AT LEAST %d GENOMES ON EACH SIDE" % MIN_PER_SIDE)
print("=" * 74)
print("  %-26s %10s %10s %10s" % ("genus", "amphibian", "endotherm", "reference"))
for g, a, o in testable:
    d = by_genus[g]
    print("  %-26s %10d %10d %10d"
          % (g, a, len(d.get("endotherm", [])), len(d.get("reference", []))))
print("  testable genera: %d" % len(testable))


def pruned(labels):
    t2 = dendropy.Tree(tree)
    t2.retain_taxa_with_labels(list(labels))
    t2.suppress_unifurcations()
    return t2


out = open(OUT, "w")
out.write("genus\tn_amphibian\tn_other\tobs_changes\tnull_mean\tnull_lo\tnull_hi\t"
          "ratio\tp\tlargest_pure_amph_clade\tverdict\n")

print()
print("=" * 74)
print("DO AMPHIBIAN AND NON-AMPHIBIAN GENOMES INTERLEAVE WITHIN EACH GENUS?")
print("=" * 74)
print("  1 change  = two clean monophyletic blocks, host nested in clade, a")
print("              within-genus gene-content test cannot separate host")
print("              from lineage.")
print("  ratio ~1  = amphibian tips scattered among the others, which is what")
print("              such a test needs.")
print()

for g, na, no in testable:
    d = by_genus[g]
    amph = d.get("amphibian", [])
    other = d.get("endotherm", []) + d.get("reference", [])
    state = {t: 1 for t in amph}
    state.update({t: 0 for t in other})
    t2 = pruned(state)
    kept = {l.taxon.label for l in t2.leaf_node_iter() if l.taxon is not None}
    state = {k: v for k, v in state.items() if k in kept}
    if len(set(state.values())) < 2:
        print("  %-26s only one state after pruning, skipped" % g)
        continue
    obs = fitch(t2, state)
    if obs < 1:
        print("  %-26s SCORE 0 WITH BOTH STATES PRESENT, traversal suspect" % g)
        continue
    labs = list(state)
    vals = [state[k] for k in labs]
    null = []
    for _ in range(N_PERM):
        random.shuffle(vals)
        null.append(fitch(t2, dict(zip(labs, vals))))
    null.sort()
    nm = sum(null) / float(len(null))
    lo = null[int(0.025 * (len(null) - 1))]
    hi = null[int(0.975 * (len(null) - 1))]
    p = (sum(1 for v in null if v <= obs) + 1) / float(N_PERM + 1)
    ratio = obs / nm if nm else float("nan")
    big = 0
    for nd in t2.postorder_internal_node_iter():
        ls = [l.taxon.label for l in nd.leaf_iter() if l.taxon is not None]
        ls = [l for l in ls if l in state]
        if len(ls) > 1 and all(state[l] == 1 for l in ls):
            big = max(big, len(ls))
    if obs <= 2:
        verdict = "separate blocks, NOT testable"
    elif ratio < 0.35:
        verdict = "strongly clustered, weak"
    elif ratio < 0.65:
        verdict = "partly interleaved"
    else:
        verdict = "interleaved, testable"
    print("  %-26s amphibian %3d vs other %3d" % (g, na, no))
    print("      changes %3d | null %6.1f [%d, %d] | ratio %.3f | p %.4f | max pure clade %d"
          % (obs, nm, lo, hi, ratio, p, big))
    print("      -> %s" % verdict)
    out.write("%s\t%d\t%d\t%d\t%.2f\t%d\t%d\t%.4f\t%.4f\t%d\t%s\n"
              % (g, na, no, obs, nm, lo, hi, ratio, p, big, verdict))

print()
print("=" * 74)
print("WHOLE TREE, AMPHIBIAN vs EVERYTHING ELSE")
print("=" * 74)
state = {}
for t in tips:
    m = meta.get(t)
    if not m or m["arm"] == "outgroup":
        continue
    state[t] = 1 if m["arm"] in AMPH_ARMS else 0
t2 = pruned(state)
kept = {l.taxon.label for l in t2.leaf_node_iter() if l.taxon is not None}
state = {k: v for k, v in state.items() if k in kept}
obs = fitch(t2, state)
labs = list(state)
vals = [state[k] for k in labs]
null = []
for _ in range(N_PERM):
    random.shuffle(vals)
    null.append(fitch(t2, dict(zip(labs, vals))))
null.sort()
nm = sum(null) / float(len(null))
p = (sum(1 for v in null if v <= obs) + 1) / float(N_PERM + 1)
na = sum(state.values())
print("  tips: %d amphibian, %d other" % (na, len(state) - na))
print("  changes %d | null mean %.1f [%d, %d] | ratio %.4f | p %.4f"
      % (obs, nm, null[0], null[-1], obs / nm if nm else float("nan"), p))
print()
print("  This is how many independent times an amphibian-associated lineage")
print("  arises. Near 1 means a single origin, no replication, and")
print("  convergence methods (treeWAS, Hogwash) have no power. Larger means")
print("  repeated origins and those methods become usable.")
out.write("WHOLE_TREE\t%d\t%d\t%d\t%.2f\t%d\t%d\t%.4f\t%.4f\t\t\n"
          % (na, len(state) - na, obs, nm, null[0], null[-1],
             obs / nm if nm else 0, p))
out.close()
print()
print("wrote %s" % OUT)

# SENTINEL_END
