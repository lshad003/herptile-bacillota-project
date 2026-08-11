# MMseqs2 clusters are merged into orthologous groups, since raw clusters split single genes.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/merge_clusters_by_og.py
# Output: work/focal_genus_pangenome/matrices/presence_og_bacteria.tsv, presence_og_narrow.tsv, cluster_to_og.tsv
# MERGE_CLUSTERS_BY_OG_V1_20260805
import os
from collections import Counter, defaultdict

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
WORK = os.path.join(ROOT, "work/focal_genus_pangenome")
MATD = os.path.join(WORK, "matrices")
CLU = os.path.join(WORK, "clu_covmode1_cluster.tsv")
ANN = os.path.join(WORK, "eggnog/focal.emapper.annotations")
MET = os.path.join(MATD, "unit_metadata.tsv")

OUT_NARROW = os.path.join(MATD, "presence_og_narrow.tsv")
OUT_BACT = os.path.join(MATD, "presence_og_bacteria.tsv")
OUT_MAP = os.path.join(MATD, "cluster_to_og.tsv")

MIN_PREV = 0.10

# Universal or near-universal genes. If merging works, these should sit near
# 1.00 in BOTH arms. In the unmerged id50_cov1 matrix rpsU was 0.76 vs 0.58
# and atpG 0.48 vs 0.19, which is impossible for genes every genome carries
# and is the signature of one gene split across several clusters.
DIAG = ["rpsU", "atpG", "rpoB", "gyrA", "recA", "dnaK", "tuf", "rpsB",
        "rplB", "frr", "pyrG", "infC", "sun", "apt", "ddh", "gpr", "cel", "prkC"]


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


meta = {r["genome"]: r for r in read_tsv(MET)}
amph = sorted(g for g, r in meta.items() if r["group"] == "amphibian")
ref = sorted(g for g, r in meta.items() if r["group"] == "reference")
allg = sorted(meta)
print("units: %d amphibian, %d reference" % (len(amph), len(ref)))

# ------------------------------------------------ cluster membership
members = defaultdict(set)
bad = 0
with open(CLU) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        g = f[1].split("|", 1)[0]
        if g in meta:
            members[f[0]].add(g)
        else:
            bad += 1
print("clusters: %d  (protein rows skipped, dropped genome: %d)" % (len(members), bad))

# ------------------------------------------------ annotations
ann = {}
hdr = None
with open(ANN) as fh:
    for line in fh:
        if line.startswith("#query"):
            hdr = line.lstrip("#").rstrip("\n").split("\t")
            continue
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if hdr and len(f) >= len(hdr):
            ann[f[0]] = dict(zip(hdr, f))
print("annotated representatives: %d" % len(ann))


def og_fields(cid):
    """eggNOG_OGs is broad-to-narrow: OG@taxid|name,OG@taxid|name,...
    Return (narrowest OG, Bacteria-level OG)."""
    a = ann.get(cid)
    if not a:
        return None, None
    s = a.get("eggNOG_OGs", "").strip()
    if not s or s == "-":
        return None, None
    parts = [p for p in s.split(",") if "@" in p]
    if not parts:
        return None, None
    narrow = parts[-1].split("@")[0]
    bact = None
    for p in parts:
        og, rest = p.split("@", 1)
        if rest.split("|")[0] == "2":
            bact = og
    if bact is None:
        bact = parts[0].split("@")[0]
    return narrow, bact


n_ann = n_noann = 0
map_narrow, map_bact = {}, {}
for cid in members:
    nog, bog = og_fields(cid)
    if nog is None:
        n_noann += 1
        # unannotated clusters keep their own identity rather than being
        # dropped, so novel gene families are not silently discarded
        map_narrow[cid] = "UNANN_" + cid
        map_bact[cid] = "UNANN_" + cid
    else:
        n_ann += 1
        map_narrow[cid] = nog
        map_bact[cid] = bog
print("clusters with an eggNOG OG: %d, without: %d" % (n_ann, n_noann))

with open(OUT_MAP, "w") as f:
    f.write("cluster\tog_narrow\tog_bacteria\tn_genomes\n")
    for cid in sorted(members):
        f.write("%s\t%s\t%s\t%d\n"
                % (cid, map_narrow[cid], map_bact[cid], len(members[cid])))
print("wrote %s" % OUT_MAP)


def build(mapping, label, outpath):
    og_members = defaultdict(set)
    og_nclust = Counter()
    for cid, gs in members.items():
        og = mapping[cid]
        og_members[og] |= gs
        og_nclust[og] += 1
    print()
    print("=" * 74)
    print("MERGED AT %s LEVEL" % label)
    print("=" * 74)
    print("  clusters %d -> orthologous groups %d" % (len(members), len(og_members)))
    multi = sum(1 for v in og_nclust.values() if v > 1)
    print("  groups built from more than one cluster: %d" % multi)
    sizes = sorted(og_nclust.values(), reverse=True)
    print("  clusters per group: max %d, median %d" % (sizes[0], sizes[len(sizes) // 2]))

    na, nr = len(amph), len(ref)
    kept = []
    for og, gs in og_members.items():
        pa = sum(1 for g in amph if g in gs) / float(na)
        pr = sum(1 for g in ref if g in gs) / float(nr)
        if max(pa, pr) >= MIN_PREV:
            kept.append(og)
    print("  groups passing the %.0f%% prevalence filter: %d" % (100 * MIN_PREV, len(kept)))

    core_a = sum(1 for og in kept
                 if sum(1 for g in amph if g in og_members[og]) / float(na) >= 0.95)
    core_r = sum(1 for og in kept
                 if sum(1 for g in ref if g in og_members[og]) / float(nr) >= 0.95)
    print("  present in >=95%% of amphibian units: %d  (was 1 unmerged)" % core_a)
    print("  present in >=95%% of reference units: %d  (was 7 unmerged)" % core_r)

    with open(outpath, "w") as f:
        f.write("cluster\t" + "\t".join(allg) + "\n")
        for og in sorted(kept):
            gs = og_members[og]
            f.write(og + "\t" + "\t".join("1" if g in gs else "0" for g in allg) + "\n")
    print("  wrote %s" % outpath)
    return og_members, mapping


def diagnostic(og_members, mapping, label):
    print()
    print("  DIAGNOSTIC: universal genes should be near 1.00 in BOTH arms")
    name_to_og = defaultdict(set)
    for cid in members:
        a = ann.get(cid)
        if not a:
            continue
        nm = a.get("Preferred_name", "").strip()
        if nm and nm != "-":
            name_to_og[nm].add(mapping[cid])
    print("    %-10s %8s %12s %12s" % ("gene", "groups", "amphibian", "reference"))
    for nm in DIAG:
        ogs = name_to_og.get(nm)
        if not ogs:
            continue
        gs = set()
        for og in ogs:
            gs |= og_members[og]
        pa = sum(1 for g in amph if g in gs) / float(len(amph))
        pr = sum(1 for g in ref if g in gs) / float(len(ref))
        print("    %-10s %8d %12.2f %12.2f" % (nm, len(ogs), pa, pr))
    print("    'groups' is how many orthologous groups carry that gene name.")
    print("    More than 1 means the merge did not fully collapse the family.")


ogm_n, _ = build(map_narrow, "NARROWEST", OUT_NARROW)
diagnostic(ogm_n, map_narrow, "narrow")

ogm_b, _ = build(map_bact, "BACTERIA", OUT_BACT)
diagnostic(ogm_b, map_bact, "bacteria")

print()
print("=" * 74)
print("HOW TO READ THIS")
print("=" * 74)
print("  If rpsU, atpG, rpoB and the ribosomal proteins now sit near 1.00 in")
print("  both arms, the merge worked and the prevalence comparison can be")
print("  rerun on the merged matrix. If they are still low, fragmentation is")
print("  worse than annotation can repair and the honest move is to report")
print("  direction only, stating that absolute prevalence is deflated.")
print("  The Bacteria level merges more aggressively and risks pooling")
print("  paralogs; the narrowest level is more conservative.")

# SENTINEL_END
