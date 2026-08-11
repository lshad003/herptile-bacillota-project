# Query alignments are assembled from the dereplication tables.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/build_rep_tree_msa.py
# Output: work/rep_tree/rep_bac120.faa, rep_tree_metadata.tsv
# BUILD_REP_TREE_MSA_V1_20260804
import os, gzip, sys
from collections import Counter, defaultdict

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
WORK = os.path.join(ROOT, "work/pooled_drep_rum")
CDB = os.path.join(WORK, "drep_out/data_tables/Cdb.csv")
WDB = os.path.join(WORK, "drep_out/data_tables/Wdb.csv")
MAP = os.path.join(WORK, "genome_arms.tsv")

OUTD = os.path.join(ROOT, "work/rep_tree")
OUT_FA = os.path.join(OUTD, "rep_bac120.faa")
OUT_META = os.path.join(OUTD, "rep_tree_metadata.tsv")

# One user MSA per arm. All produced by GTDB-Tk against r220 with the
# canonical mask, so columns are positionally homologous.
MSAS = {
    "herptile": "results/gtdbtk_wild_sgb_r220/align/gtdbtk.bac120.user_msa.fasta.gz",
    "youngblut": "results/gtdbtk_youngblut_r220/align/gtdbtk.bac120.user_msa.fasta.gz",
    "ehi": "results/gtdbtk_ehi_r220/align/gtdbtk.bac120.user_msa.fasta.gz",
    "ehi_amphibian": "results/gtdbtk_ehi_amphibian_r220/align/gtdbtk.bac120.user_msa.fasta.gz",
}
SUMS = {
    "herptile": "results/gtdbtk_wild_sgb_r220/gtdbtk.bac120.summary.tsv",
    "youngblut": "results/gtdbtk_youngblut_r220/gtdbtk.bac120.summary.tsv",
    "ehi": "results/gtdbtk_ehi_r220_classify/gtdbtk.bac120.summary.tsv",
    "ehi_amphibian": "results/gtdbtk_ehi_amphibian_r220/gtdbtk.bac120.summary.tsv",
}

SGB_MANIFEST = os.path.join(ROOT, "data/sgb_manifest.tsv")
HERP_MANIFEST = os.path.join(ROOT, "data/herptile_bacillota_A_HQ_manifest_with_source.tsv")
EHI_MANIFEST = os.path.join(ROOT, "results/ehi_nonherptile_manifest.tsv")
EHI_AMPH_MANIFEST = os.path.join(ROOT, "results/ehi_amphibian_manifest.tsv")
YB_QC = os.path.join(ROOT, "data/youngblut/youngblut_fetch_qc.tsv")
GTDB_TAX = "/srv/projects/db/gtdbtk/220/taxonomy/gtdb_taxonomy.tsv"

# gtdb_ref genomes have no user MSA: they are already IN the r220 reference
# alignment. Including them would need the reference MSA, which is a
# different and much larger file. They are excluded from this tree and that
# is a stated limitation, not an oversight.
SKIP_ARMS = {"gtdb_ref"}

AMPHIBIAN = {"Salamander", "Frog", "Toad"}
REPTILE = {"Lizard", "Turtle", "Tortoise", "Snake"}


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


def read_fasta_gz(path):
    seqs = {}
    name, buf = None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


if not os.path.isdir(OUTD):
    os.makedirs(OUTD)

# ------------------------------------------------------------ arms + reps
arm_of, gid_of = {}, {}
with open(MAP) as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        arm_of[f[0]] = f[1]
        gid_of[f[0]] = f[2]
print("=" * 74)
print("POOLED dRep")
print("=" * 74)
print("  genomes mapped: %d  %s" % (len(arm_of), dict(Counter(arm_of.values()))))

clusters = defaultdict(list)
with open(CDB) as fh:
    hdr = fh.readline().rstrip("\n").replace('"', "").split(",")
    gi, si = hdr.index("genome"), hdr.index("secondary_cluster")
    for line in fh:
        f = line.rstrip("\n").replace('"', "").split(",")
        if len(f) > max(gi, si):
            clusters[f[si]].append(f[gi])
print("  clusters: %d" % len(clusters))

# ONE TIP PER ARM PER CLUSTER, not one per cluster. A cluster containing
# genomes from two arms is informative and must not be collapsed to one.
tips = []
for cid, gs in clusters.items():
    per_arm = defaultdict(list)
    for g in gs:
        per_arm[arm_of.get(g, "?")].append(g)
    for a, members in per_arm.items():
        if a in SKIP_ARMS:
            continue
        # pick the member with the most alignment later; for now take the
        # lexicographically first so the choice is deterministic
        tips.append((cid, a, sorted(members)[0], len(members)))
print("  tips (one per arm per cluster, gtdb_ref excluded): %d" % len(tips))
print("  by arm: %s" % dict(Counter(a for _, a, _, _ in tips)))

# ------------------------------------------------------------ load MSAs
print()
print("=" * 74)
print("USER MSAs")
print("=" * 74)
msa, widths = {}, {}
for arm, rel in MSAS.items():
    p = os.path.join(ROOT, rel)
    if not os.path.exists(p):
        print("  %-14s MISSING (%s)" % (arm, rel))
        continue
    s = read_fasta_gz(p)
    w = set(len(v) for v in s.values())
    msa[arm] = s
    widths[arm] = w
    print("  %-14s %5d seqs, widths %s" % (arm, len(s), w))

if not msa:
    sys.exit("no MSAs found")
allw = set()
for w in widths.values():
    allw |= w
if len(allw) != 1:
    sys.exit("MSA widths differ across arms: %s. Cannot concatenate." % allw)
W = allw.pop()
print("  shared width %d, concatenation is column-consistent" % W)

missing_arms = [a for _, a, _, _ in tips if a not in msa]
if missing_arms:
    print()
    print("  ARMS WITH TIPS BUT NO MSA: %s" % dict(Counter(missing_arms)))
    print("  Those tips will be dropped. If ehi_amphibian is listed, the")
    print("  classify job has not finished; rerun this after it lands.")

# ------------------------------------------------------------ metadata
tax = {}
for arm, rel in SUMS.items():
    p = os.path.join(ROOT, rel)
    if not os.path.exists(p):
        continue
    for r in read_tsv(p):
        tax[(arm, r["user_genome"].strip())] = parse_tax(r["classification"])

taxon_map = {}
for r in read_tsv(HERP_MANIFEST):
    t = r["host_taxon"].strip()
    if t:
        taxon_map[t] = r["animal_type"].strip()

herp_host = {}
for r in read_tsv(SGB_MANIFEST):
    if r["has_wild"].strip().lower() != "yes":
        continue
    g = {taxon_map[t] for t in r["host_species"].split(";")
         if t.strip() and t.strip() in taxon_map}
    herp_host[r["representative"].strip()] = g

ehi_host, ehi_class = {}, {}
for r in read_tsv(EHI_MANIFEST):
    ehi_host[r["genome_id"].strip()] = r["host_species"].strip()
    ehi_class[r["genome_id"].strip()] = r["host_class"].strip()

amph_host = {}
if os.path.exists(EHI_AMPH_MANIFEST):
    for r in read_tsv(EHI_AMPH_MANIFEST):
        amph_host[r["genome_id"].strip()] = r["host_species"].strip()

yb_host = {}
for r in read_tsv(YB_QC):
    yb_host[r["genome"].strip()] = r["host"].strip()


def host_info(arm, gid):
    if arm == "herptile":
        g = herp_host.get(gid, set())
        if g & AMPHIBIAN and g & REPTILE:
            cat = "mixed"
        elif g & AMPHIBIAN:
            cat = "amphibian"
        elif g & REPTILE:
            cat = "reptile"
        else:
            cat = "unknown"
        return cat, ";".join(sorted(g))
    if arm == "ehi":
        return "endotherm", ehi_host.get(gid, "")
    if arm == "ehi_amphibian":
        return "amphibian", amph_host.get(gid, "")
    if arm == "youngblut":
        return "endotherm", yb_host.get(gid, "")
    return "unknown", ""


# ------------------------------------------------------------ write
written = 0
skipped = 0
with open(OUT_FA, "w") as fa, open(OUT_META, "w") as mp:
    mp.write("tip\tcluster\tarm\tgenome\thost_category\thost_detail\t"
             "family\tgenus\tn_in_cluster_this_arm\tgap_fraction\n")
    for cid, a, g, n in sorted(tips):
        gid = gid_of.get(g, "")
        s = msa.get(a, {}).get(gid)
        if s is None:
            skipped += 1
            continue
        tip = "%s|%s|%s" % (a, cid, gid)
        cat, detail = host_info(a, gid)
        t = tax.get((a, gid), {})
        gf = s.count("-") / float(len(s))
        fa.write(">%s\n%s\n" % (tip, s))
        mp.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%.4f\n"
                 % (tip, cid, a, gid, cat, detail,
                    t.get("f", ""), t.get("g", ""), n, gf))
        written += 1

print()
print("=" * 74)
print("OUTPUT")
print("=" * 74)
print("  tips written: %d, skipped for missing sequence: %d" % (written, skipped))
print("  wrote %s" % OUT_FA)
print("  wrote %s" % OUT_META)

rows = read_tsv(OUT_META)
print()
print("  by arm          : %s" % dict(Counter(r["arm"] for r in rows)))
print("  by host category: %s" % dict(Counter(r["host_category"] for r in rows)))
gf = defaultdict(list)
for r in rows:
    gf[r["arm"]].append(float(r["gap_fraction"]))
print()
print("  gap fraction by arm (arm-correlated gappiness drives long-branch")
print("  artefacts, so this must be checked before any inference):")
for a in sorted(gf):
    v = sorted(gf[a])
    print("    %-14s n=%4d  min %.3f median %.3f max %.3f"
          % (a, len(v), v[0], v[len(v) // 2], v[-1]))

mixed = 0
for cid, gs in clusters.items():
    arms = set(arm_of.get(g, "?") for g in gs) - SKIP_ARMS
    if len(arms) > 1:
        mixed += 1
print()
print("  clusters spanning more than one non-reference arm: %d" % mixed)
print("  Those are the informative ones: a cluster holding both an amphibian")
print("  and an endotherm genome breaks the host/clade confound locally, which")
print("  is what the gate diagnostic needs to have any power.")

# SENTINEL_END
