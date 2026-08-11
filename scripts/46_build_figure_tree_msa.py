# Query and reference alignments are combined with an outgroup into the tree input.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/build_figure_tree_msa.py
# Output: work/rep_tree/figure_tree.faa, figure_tree_metadata.tsv
# BUILD_FIGURE_TREE_MSA_V1_20260804
import os, gzip, sys
from collections import Counter, defaultdict

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
WORK = os.path.join(ROOT, "work/pooled_drep_rum")
CDB = os.path.join(WORK, "drep_out/data_tables/Cdb.csv")
MAP = os.path.join(WORK, "genome_arms.tsv")
OUTD = os.path.join(ROOT, "work/rep_tree")
OUT_FA = os.path.join(OUTD, "figure_tree.faa")
OUT_META = os.path.join(OUTD, "figure_tree_metadata.tsv")

REF_MSA = "/srv/projects/db/gtdbtk/220/msa/gtdb_r220_bac120.faa"
MASK = "/srv/projects/db/gtdbtk/220/masks/gtdb_r220_bac120.mask"
GTDB_TAX = "/srv/projects/db/gtdbtk/220/taxonomy/gtdb_taxonomy.tsv"
REF_MASKED = os.path.join(OUTD, "gtdb_ref_rum_masked.faa")

MSAS = {
    "herptile": "results/gtdbtk_wild_sgb_r220/align/gtdbtk.bac120.user_msa.fasta.gz",
    "ehi_amphibian": "results/gtdbtk_ehi_amphibian_r220/align/gtdbtk.bac120.user_msa.fasta.gz",
    "ehi": "results/gtdbtk_ehi_r220/align/gtdbtk.bac120.user_msa.fasta.gz",
    "youngblut": "results/gtdbtk_youngblut_r220/align/gtdbtk.bac120.user_msa.fasta.gz",
}
SUMS = {
    "herptile": "results/gtdbtk_wild_sgb_r220/gtdbtk.bac120.summary.tsv",
    "ehi_amphibian": "results/gtdbtk_ehi_amphibian_r220/gtdbtk.bac120.summary.tsv",
    "ehi": "results/gtdbtk_ehi_r220_classify/gtdbtk.bac120.summary.tsv",
    "youngblut": "results/gtdbtk_youngblut_r220/gtdbtk.bac120.summary.tsv",
}
SGB_MANIFEST = os.path.join(ROOT, "data/sgb_manifest.tsv")
HERP_MANIFEST = os.path.join(ROOT, "data/herptile_bacillota_A_HQ_manifest_with_source.tsv")
AMPH_MANIFEST = os.path.join(ROOT, "results/ehi_amphibian_manifest.tsv")
EHI_MANIFEST = os.path.join(ROOT, "results/ehi_nonherptile_manifest.tsv")
YB_QC = os.path.join(ROOT, "data/youngblut/youngblut_fetch_qc.tsv")

OUTGROUP_FAMS = ("Acutalibacteraceae", "Butyricicoccaceae")
N_OUTGROUP_PER_FAM = 6
W = 5035
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


def read_fa(path, opener=open):
    seqs, name, buf = {}, None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


if not os.path.isdir(OUTD):
    os.makedirs(OUTD)

arm_of, gid_of = {}, {}
with open(MAP) as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        arm_of[f[0]] = f[1]
        gid_of[f[0]] = f[2]

clusters = defaultdict(list)
with open(CDB) as fh:
    hdr = fh.readline().rstrip("\n").replace('"', "").split(",")
    gi, si = hdr.index("genome"), hdr.index("secondary_cluster")
    for line in fh:
        f = line.rstrip("\n").replace('"', "").split(",")
        if len(f) > max(gi, si):
            clusters[f[si]].append(f[gi])
print("=" * 74)
print("FIVE-ARM dRep")
print("=" * 74)
print("  clusters: %d, genomes: %d"
      % (len(clusters), sum(len(v) for v in clusters.values())))

tips = []
for cid, gs in clusters.items():
    per_arm = defaultdict(list)
    for g in gs:
        per_arm[arm_of.get(g, "?")].append(g)
    for a, members in per_arm.items():
        if a == "gtdb_ref":
            continue
        tips.append((cid, a, sorted(members)[0]))
print("  query tips (one per arm per cluster): %d" % len(tips))
print("  by arm: %s" % dict(Counter(a for _, a, _ in tips)))

msa = {}
for arm, rel in MSAS.items():
    p = os.path.join(ROOT, rel)
    if not os.path.exists(p):
        sys.exit("missing MSA: %s" % p)
    msa[arm] = read_fa(p, gzip.open)
    w = set(len(v) for v in msa[arm].values())
    print("  %-14s %5d seqs, widths %s" % (arm, len(msa[arm]), w))
    if w != {W}:
        sys.exit("width mismatch in %s" % arm)

if not os.path.exists(REF_MASKED):
    sys.exit("run extract_ref_msa_rum.py first: %s missing" % REF_MASKED)
refs = read_fa(REF_MASKED)
print("  gtdb refs (masked): %d" % len(refs))

print()
print("=" * 74)
print("OUTGROUP: %s" % ", ".join(OUTGROUP_FAMS))
print("=" * 74)
og_want = defaultdict(list)
with open(GTDB_TAX) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        fam = parse_tax(f[1]).get("f", "")
        if fam in OUTGROUP_FAMS:
            og_want[fam].append(f[0].strip().replace("GB_", "").replace("RS_", ""))
for fam in OUTGROUP_FAMS:
    print("  %-24s in GTDB: %d" % (fam, len(og_want[fam])))
og_target = set()
for fam in OUTGROUP_FAMS:
    og_target |= set(sorted(og_want[fam])[:N_OUTGROUP_PER_FAM])
print("  selected: %d" % len(og_target))
if not og_target:
    sys.exit("OUTGROUP EMPTY: family names do not match GTDB. Tree would be unrooted.")

with open(MASK) as fh:
    mask = fh.readline().strip()
keep = [i for i, c in enumerate(mask) if c == "1"]
if len(keep) != W:
    sys.exit("mask has %d ones, expected %d" % (len(keep), W))

og = {}
cur, buf, on = None, [], False
with open(REF_MSA) as fh:
    for line in fh:
        if line.startswith(">"):
            if cur and on:
                s = "".join(buf)
                if len(s) == len(mask):
                    og[cur] = "".join(s[i] for i in keep)
            cur = line[1:].strip().split()[0].replace("GB_", "").replace("RS_", "")
            on = cur in og_target
            buf = []
        elif on:
            buf.append(line.strip())
    if cur and on:
        s = "".join(buf)
        if len(s) == len(mask):
            og[cur] = "".join(s[i] for i in keep)
print("  outgroup sequences masked: %d" % len(og))

tax = {}
for arm, rel in SUMS.items():
    for r in read_tsv(os.path.join(ROOT, rel)):
        tax[(arm, r["user_genome"].strip())] = parse_tax(r["classification"])

taxon_map = {}
for r in read_tsv(HERP_MANIFEST):
    if r["host_taxon"].strip():
        taxon_map[r["host_taxon"].strip()] = r["animal_type"].strip()
herp_host = {}
for r in read_tsv(SGB_MANIFEST):
    if r["has_wild"].strip().lower() == "yes":
        herp_host[r["representative"].strip()] = {
            taxon_map[t] for t in r["host_species"].split(";")
            if t.strip() and t.strip() in taxon_map}
amph_host = {r["genome_id"].strip(): r["host_species"].strip()
             for r in read_tsv(AMPH_MANIFEST)}
ehi_host = {r["genome_id"].strip(): r["host_species"].strip()
            for r in read_tsv(EHI_MANIFEST)}
yb_host = {r["genome"].strip(): r["host"].strip() for r in read_tsv(YB_QC)}

ref_genus = {}
with open(GTDB_TAX) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        t = parse_tax(f[1])
        if t.get("f", "") == "Ruminococcaceae":
            ref_genus[f[0].strip().replace("GB_", "").replace("RS_", "")] = t.get("g", "")


def host_of(arm, gid):
    if arm == "herptile":
        g = herp_host.get(gid, set())
        if g & AMPHIBIAN and not g & REPTILE:
            cat = "amphibian"
        elif g & REPTILE and not g & AMPHIBIAN:
            cat = "reptile"
        elif g:
            cat = "mixed"
        else:
            cat = "unknown"
        return cat, ";".join(sorted(g))
    if arm == "ehi_amphibian":
        return "amphibian", amph_host.get(gid, "")
    if arm == "ehi":
        return "endotherm", ehi_host.get(gid, "")
    if arm == "youngblut":
        return "endotherm", yb_host.get(gid, "")
    return "none", ""


n = 0
with open(OUT_FA, "w") as fa, open(OUT_META, "w") as mp:
    mp.write("tip\tarm\tgenome\tcluster\thost_category\thost_detail\t"
             "family\tgenus\tgap_fraction\n")
    for cid, a, g in sorted(tips):
        gid = gid_of.get(g, "")
        s = msa.get(a, {}).get(gid)
        if s is None:
            continue
        cat, detail = host_of(a, gid)
        t = tax.get((a, gid), {})
        tip = "%s|%s" % (a, gid)
        fa.write(">%s\n%s\n" % (tip, s))
        mp.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.4f\n"
                 % (tip, a, gid, cid, cat, detail, t.get("f", ""),
                    t.get("g", ""), s.count("-") / float(W)))
        n += 1
    for acc, s in sorted(refs.items()):
        tip = "gtdb_ref|%s" % acc
        fa.write(">%s\n%s\n" % (tip, s))
        mp.write("%s\tgtdb_ref\t%s\t\tnone\t\tRuminococcaceae\t%s\t%.4f\n"
                 % (tip, acc, ref_genus.get(acc, ""), s.count("-") / float(W)))
        n += 1
    for acc, s in sorted(og.items()):
        tip = "outgroup|%s" % acc
        fa.write(">%s\n%s\n" % (tip, s))
        mp.write("%s\toutgroup\t%s\t\tnone\t\toutgroup\t\t%.4f\n"
                 % (tip, acc, s.count("-") / float(W)))
        n += 1

print()
print("=" * 74)
print("OUTPUT")
print("=" * 74)
print("  total tips: %d" % n)
rows = read_tsv(OUT_META)
print("  by arm: %s" % dict(Counter(r["arm"] for r in rows)))
print("  by host: %s" % dict(Counter(r["host_category"] for r in rows)))
gf = defaultdict(list)
for r in rows:
    gf[r["arm"]].append(float(r["gap_fraction"]))
print()
print("  gap fraction by arm:")
for a in sorted(gf):
    v = sorted(gf[a])
    print("    %-14s n=%4d min %.3f median %.3f max %.3f"
          % (a, len(v), v[0], v[len(v) // 2], v[-1]))
print()
print("  wrote %s" % OUT_FA)
print("  wrote %s" % OUT_META)

# SENTINEL_END
