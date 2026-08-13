# Genus-unassigned query tips are grouped into monophyletic clades and their sister genera identified.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/unassigned_clade_coherence.py
# Output: results/unassigned_clade_coherence.tsv
import os, sys, csv
from collections import Counter, defaultdict
import dendropy

BASE = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
TREE = os.path.join(BASE, "work/rep_tree/figure_tree.nwk")
META = os.path.join(BASE, "work/rep_tree/figure_tree_metadata_genus.tsv")
AUDIT = os.path.join(BASE, "results/seqcode_representative_audit_v2.tsv")
OUT  = os.path.join(BASE, "results/unassigned_clade_coherence.tsv")

QUERY_ARMS = ("herptile", "ehi_amphibian")
UNASSIGNED = ("", "UNASSIGNED", "unassigned", "g__", "NA", "na", "None", "-")

def die(m):
    print("")
    print("!" * 72)
    print("FAILED: " + m)
    print("!" * 72)
    sys.exit(1)

def sniff(p):
    with open(p) as fh:
        return "," if fh.readline().count(",") > 1 else "\t"

def norm(s):
    s = str(s).strip()
    if "|" in s:
        s = s.split("|")[-1]
    return s

print("=" * 72)
print("STEP 1  TREE AND METADATA")
print("=" * 72)
for p in (TREE, META):
    if not os.path.exists(p):
        die("missing " + p)
mrows = list(csv.DictReader(open(META), delimiter=sniff(META)))
meta = {}
for r in mrows:
    meta[r["tip"].strip()] = {"arm": r["arm"].strip(),
                              "genus": (r["genus"] or "").strip(),
                              "genome": norm(r["genome"]),
                              "family": (r.get("family") or "").strip()}
tree = dendropy.Tree.get(path=TREE, schema="newick", preserve_underscores=True)
leaves = list(tree.leaf_node_iter())
labels = [nd.taxon.label if nd.taxon is not None else "" for nd in leaves]
n = len(leaves)
print("  tips: %d   metadata rows: %d" % (n, len(mrows)))
if sum(1 for L in labels if L in meta) != n:
    die("metadata does not cover every tip")
print("  arms: %s" % dict(Counter(meta[L]["arm"] for L in labels)))

og = [nd for nd in leaves if meta[nd.taxon.label]["arm"] == "outgroup"]
print("  outgroup tips: %d" % len(og))
if len(og) >= 2:
    mr = tree.mrca(taxon_labels=[nd.taxon.label for nd in og])
    under = len(list(mr.leaf_iter()))
    print("  outgroup MRCA subtends %d tips (%d expected if monophyletic)" % (under, len(og)))
    if under == len(og):
        tree.reroot_at_edge(mr.edge, update_bipartitions=False)
        print("  REROOTED on the outgroup")
    else:
        print("  outgroup NOT monophyletic, tree left as read. Monophyly calls below")
        print("  are therefore on the tree as given, state that in any caption.")
tree.calc_node_root_distances(return_leaf_distances_only=False)
rd = {nd: (nd.root_distance or 0.0) for nd in tree.preorder_node_iter()}
leaves = list(tree.leaf_node_iter())
labels = [nd.taxon.label if nd.taxon is not None else "" for nd in leaves]
bylabel = {nd.taxon.label: nd for nd in leaves if nd.taxon is not None}

def anc_set(nd):
    s, x = set(), nd
    while x is not None:
        s.add(id(x))
        x = x.parent_node
    return s

def mrca2(a, b, acache=None):
    s = acache if acache is not None else anc_set(a)
    x = b
    while x is not None:
        if id(x) in s:
            return x
        x = x.parent_node
    return None

def pd(a, b, acache=None):
    m = mrca2(a, b, acache)
    if m is None:
        return None
    return rd[a] + rd[b] - 2 * rd[m]

def is_un(g):
    return g.strip() in UNASSIGNED

print("")
print("=" * 72)
print("STEP 2  DOES THIS TREE RECOVER NAMED GENERA? (the calibration)")
print("=" * 72)
bygen = defaultdict(list)
for L in labels:
    g = meta[L]["genus"]
    if not is_un(g):
        bygen[g].append(L)
multi = {g: v for g, v in bygen.items() if len(v) >= 3}
print("  named genera with >=3 tips: %d" % len(multi))
mono = nonmono = 0
nonmono_ex = []
for g, tips in multi.items():
    mr = tree.mrca(taxon_labels=tips)
    under = set(x.taxon.label for x in mr.leaf_iter() if x.taxon is not None)
    if under == set(tips):
        mono += 1
    else:
        nonmono += 1
        if len(nonmono_ex) < 6:
            intr = under - set(tips)
            nonmono_ex.append((g, len(tips), len(under), len(intr)))
print("  MONOPHYLETIC     %d (%.1f%%)" % (mono, 100.0 * mono / max(len(multi), 1)))
print("  NOT monophyletic %d (%.1f%%)" % (nonmono, 100.0 * nonmono / max(len(multi), 1)))
for g, nt, nu, ni in nonmono_ex:
    print("     %-24s %d tips, MRCA subtends %d, %d intruders" % (g, nt, nu, ni))
print("  READ THIS: if named genera are largely monophyletic here, a monophyly")
print("  call on an unassigned group is meaningful. If they are not, the tree")
print("  cannot support genus delineation and nothing below should be used.")

print("")
print("=" * 72)
print("STEP 3  UNASSIGNED QUERY TIPS")
print("=" * 72)
un_tips = [L for L in labels
           if meta[L]["arm"] in QUERY_ARMS and is_un(meta[L]["genus"])]
print("  query tips total    : %d" % sum(1 for L in labels if meta[L]["arm"] in QUERY_ARMS))
print("  of which UNASSIGNED : %d" % len(un_tips))
print("  by arm: %s" % dict(Counter(meta[L]["arm"] for L in un_tips)))
if not un_tips:
    die("no unassigned query tips on this tree, nothing to test")

unset = set(un_tips)
qual = []
for nd in tree.preorder_node_iter():
    lv = [x.taxon.label for x in nd.leaf_iter() if x.taxon is not None]
    if lv and all(L in unset for L in lv):
        p = nd.parent_node
        pok = False
        if p is not None:
            plv = [x.taxon.label for x in p.leaf_iter() if x.taxon is not None]
            pok = all(L in unset for L in plv)
        if not pok:
            qual.append((nd, lv))
sizes = Counter(len(lv) for _, lv in qual)
print("")
print("  MAXIMAL all-unassigned clades: %d" % len(qual))
print("  size distribution: %s" % dict(sorted(sizes.items())))
mult = [(nd, lv) for nd, lv in qual if len(lv) >= 2]
sing = [(nd, lv) for nd, lv in qual if len(lv) == 1]
print("  clades with >=2 tips: %d  (covering %d tips)"
      % (len(mult), sum(len(lv) for _, lv in mult)))
print("  singletons          : %d" % len(sing))

print("")
print("=" * 72)
print("STEP 4  WITHIN-GENUS DISTANCE, EMPIRICAL CALIBRATION")
print("=" * 72)
wg = []
for g, tips in multi.items():
    if len(tips) > 25:
        tips = tips[:25]
    ds = []
    for i in range(len(tips)):
        a = bylabel[tips[i]]
        ac = anc_set(a)
        for j in range(i + 1, len(tips)):
            d = pd(a, bylabel[tips[j]], ac)
            if d is not None:
                ds.append(d)
    if ds:
        wg.append((g, sum(ds) / len(ds), max(ds)))
if wg:
    means = sorted(x[1] for x in wg)
    maxes = sorted(x[2] for x in wg)
    def q(v, f):
        return v[min(len(v) - 1, int(f * len(v)))]
    print("  named genera measured: %d" % len(wg))
    print("  MEAN within-genus patristic distance:")
    print("    median %.4f   10th %.4f   90th %.4f" % (q(means, .5), q(means, .1), q(means, .9)))
    print("  MAX within-genus patristic distance:")
    print("    median %.4f   90th %.4f   max %.4f" % (q(maxes, .5), q(maxes, .9), maxes[-1]))
    print("  THIS IS DESCRIPTIVE. It is NOT a cutoff and must not become one.")
    print("  Genus boundaries here come from GTDB-Tk against the full database.")

print("")
print("=" * 72)
print("STEP 5  EACH CANDIDATE CLADE, WITH ITS SISTER")
print("=" * 72)
pass4 = set()
if os.path.exists(AUDIT):
    for r in csv.DictReader(open(AUDIT), delimiter=sniff(AUDIT)):
        if r.get("pass_four", "0").strip() == "1":
            pass4.add(r["representative"].strip())
    print("  representatives passing all four assembly criteria: %d" % len(pass4))
else:
    print("  audit file absent, pass_four column will read NA")

rows = []
mult.sort(key=lambda t: -len(t[1]))
for k, (nd, lv) in enumerate(mult, 1):
    p = nd.parent_node
    sis = []
    if p is not None:
        for c in p.child_node_iter():
            if c is nd:
                continue
            sis += [x.taxon.label for x in c.leaf_iter() if x.taxon is not None]
    sg = Counter(meta[L]["genus"] if not is_un(meta[L]["genus"]) else "UNASSIGNED" for L in sis)
    sa = Counter(meta[L]["arm"] for L in sis)
    a0 = bylabel[lv[0]]
    ac = anc_set(a0)
    intern = [pd(a0, bylabel[x], ac) for x in lv[1:]]
    intern = [d for d in intern if d is not None]
    near, nearlab = None, ""
    for L in labels:
        if L in unset or meta[L]["arm"] == "outgroup":
            continue
        d = pd(a0, bylabel[L], ac)
        if d is not None and (near is None or d < near):
            near, nearlab = d, L
    genomes = [meta[L]["genome"] for L in lv]
    npass = sum(1 for g in genomes if g in pass4)
    print("")
    print("  CLADE %d: %d tips" % (k, len(lv)))
    print("    arms      : %s" % dict(Counter(meta[L]["arm"] for L in lv)))
    print("    max internal distance: %s"
          % ("%.4f" % max(intern) if intern else "n/a"))
    print("    sister    : %d tips, arms %s" % (len(sis), dict(sa)))
    print("    sister genera: %s" % dict(sg.most_common(4)))
    if near is not None:
        print("    nearest assigned tip: %.4f (%s, genus %s)"
              % (near, nearlab, meta[nearlab]["genus"] or "UNASSIGNED"))
    print("    genomes passing all four assembly criteria: %d of %d" % (npass, len(lv)))
    for L in lv[:8]:
        g = meta[L]["genome"]
        print("       %-42s %-14s %s" % (g, meta[L]["arm"], "PASS4" if g in pass4 else ""))
    if len(lv) > 8:
        print("       ... and %d more" % (len(lv) - 8))
    rows.append((k, len(lv), max(intern) if intern else "", near if near else "",
                 nearlab, npass, ";".join(genomes)))

sp4 = sum(1 for nd, lv in sing if meta[lv[0]]["genome"] in pass4)
print("")
print("  SINGLETONS: %d, of which %d pass all four assembly criteria" % (len(sing), sp4))
print("  A singleton cannot support a genus proposal on monophyly alone.")

if os.path.exists(OUT):
    print("")
    print("  NOT overwriting existing " + OUT)
else:
    with open(OUT, "w") as fh:
        fh.write("clade\tn_tips\tmax_internal_dist\tdist_to_nearest_assigned\tnearest_tip\tn_pass_four\tgenomes\n")
        for r in rows:
            fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % r)
        for nd, lv in sing:
            g = meta[lv[0]]["genome"]
            fh.write("singleton\t1\t\t\t\t%d\t%s\n" % (1 if g in pass4 else 0, g))
    print("")
    print("  wrote " + OUT)

print("")
print("=" * 72)
print("WHAT THIS DOES AND DOES NOT ESTABLISH")
print("=" * 72)
print("  DOES: whether GTDB-Tk genus-unassigned query genomes group together on")
print("    a tree built with COMPLETE reference membership, and how many")
print("    independent groups there are.")
print("  DOES NOT: assign a rank. Genus delineation in GTDB uses RED, which")
print("    needs the calibrated reference tree and is NOT computed here.")
print("  NO PATRISTIC CUTOFF WAS APPLIED. Groups come from GTDB-Tk assignment")
print("    plus monophyly. Distances above are descriptive calibration only.")
print("    This is deliberate: a self-chosen cutoff on a filtered reference set")
print("    is what invalidated the earlier novel-clade claim.")
print("  SCOPE: this tree is Ruminococcaceae only. Unassigned representatives in")
print("    other families are NOT covered and need their own tree.")
print("")
print("UNASSIGNED_CLADE_COHERENCE_V1_20260806 COMPLETE")
# UNASSIGNED_CLADE_COHERENCE_V1_20260806
