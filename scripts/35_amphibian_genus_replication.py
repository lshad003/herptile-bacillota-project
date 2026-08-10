# Genus recovery is cross-tabulated between the two amphibian catalogs and the endotherm catalogs.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/amphibian_genus_replication.py
# Output: results/amphibian_genus_replication.tsv
# AMPHIBIAN_GENUS_REPLICATION_V1_20260804
import os
from collections import Counter, defaultdict

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"

SGB_MANIFEST = os.path.join(ROOT, "data/sgb_manifest.tsv")
HERP_MANIFEST = os.path.join(ROOT, "data/herptile_bacillota_A_HQ_manifest_with_source.tsv")
HERP_SUM = os.path.join(ROOT, "results/gtdbtk_wild_sgb_r220/gtdbtk.bac120.summary.tsv")
AMPH_SUM = os.path.join(ROOT, "results/gtdbtk_ehi_amphibian_r220/gtdbtk.bac120.summary.tsv")
AMPH_MAN = os.path.join(ROOT, "results/ehi_amphibian_manifest.tsv")
EHI_SUM = os.path.join(ROOT, "results/gtdbtk_ehi_r220_classify/gtdbtk.bac120.summary.tsv")
EHI_MAN = os.path.join(ROOT, "results/ehi_nonherptile_manifest.tsv")
YB_SUM = os.path.join(ROOT, "results/gtdbtk_youngblut_r220/gtdbtk.bac120.summary.tsv")
GTDB_TAX = "/srv/projects/db/gtdbtk/220/taxonomy/gtdb_taxonomy.tsv"

OUT = os.path.join(ROOT, "results/amphibian_genus_replication.tsv")

FAM = "Ruminococcaceae"
AMPHIBIAN = {"Salamander", "Frog", "Toad"}


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


def norm(g):
    g = g.strip()
    return "" if g.upper() in ("UNASSIGNED", "", "NA", "N/A") else g


# --------------------------------------------------- herptile, by host group
taxon_map = {}
for r in read_tsv(HERP_MANIFEST):
    t = r["host_taxon"].strip()
    if t:
        taxon_map[t] = r["animal_type"].strip()

herp_gen = {}
for r in read_tsv(HERP_SUM):
    t = parse_tax(r["classification"])
    if t.get("f", "") == FAM:
        herp_gen[r["user_genome"].strip()] = norm(t.get("g", ""))

herp_by_group = defaultdict(Counter)
herp_all = Counter()
for r in read_tsv(SGB_MANIFEST):
    if r["has_wild"].strip().lower() != "yes":
        continue
    rep = r["representative"].strip()
    g = herp_gen.get(rep)
    if not g:
        continue
    herp_all[g] += 1
    for t in r["host_species"].split(";"):
        t = t.strip()
        if t in taxon_map:
            herp_by_group[taxon_map[t]][g] += 1

# --------------------------------------------------- ehi newts, by host
amph_host = {r["genome_id"].strip(): r["host_species"].strip()
             for r in read_tsv(AMPH_MAN)}
newt_all = Counter()
newt_by_host = defaultdict(Counter)
for r in read_tsv(AMPH_SUM):
    t = parse_tax(r["classification"])
    if t.get("f", "") != FAM:
        continue
    g = norm(t.get("g", ""))
    if not g:
        continue
    gid = r["user_genome"].strip()
    newt_all[g] += 1
    newt_by_host[amph_host.get(gid, "?")][g] += 1

# --------------------------------------------------- endotherm arms
ehi_cls = {r["genome_id"].strip(): r["host_class"].strip()
           for r in read_tsv(EHI_MAN)}
ehi_mam = Counter()
for r in read_tsv(EHI_SUM):
    t = parse_tax(r["classification"])
    if t.get("f", "") != FAM:
        continue
    g = norm(t.get("g", ""))
    if g and ehi_cls.get(r["user_genome"].strip(), "") == "Mammalia":
        ehi_mam[g] += 1

yb = Counter()
for r in read_tsv(YB_SUM):
    t = parse_tax(r["classification"])
    if t.get("f", "") == FAM:
        g = norm(t.get("g", ""))
        if g:
            yb[g] += 1

ref = Counter()
with open(GTDB_TAX) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        t = parse_tax(f[1])
        if t.get("f", "") == FAM:
            g = norm(t.get("g", ""))
            if g:
                ref[g] += 1

print("=" * 74)
print("ARMS, GENUS LEVEL, %s" % FAM)
print("=" * 74)
print("  herptile wild : %4d genomes, %3d genera" % (sum(herp_all.values()), len(herp_all)))
print("  EHI newts     : %4d genomes, %3d genera" % (sum(newt_all.values()), len(newt_all)))
print("  EHI mammals   : %4d genomes, %3d genera" % (sum(ehi_mam.values()), len(ehi_mam)))
print("  Youngblut     : %4d genomes, %3d genera" % (sum(yb.values()), len(yb)))
print("  GTDB refs     : %4d genomes, %3d genera" % (sum(ref.values()), len(ref)))
print()
print("  EHI newt hosts: %s" % {k: sum(v.values()) for k, v in newt_by_host.items()})

H, N, M, Y = set(herp_all), set(newt_all), set(ehi_mam), set(yb)
endo = M | Y

print()
print("=" * 74)
print("DOES THE NEWT CATALOG REPLICATE THE HERPTILE GENUS SET?")
print("=" * 74)
print("  herptile genera            : %d" % len(H))
print("  newt genera                : %d" % len(N))
print("  shared herptile and newt   : %d" % len(H & N))
print("  jaccard herptile vs newt   : %.3f" % (len(H & N) / float(len(H | N))))
print("  jaccard herptile vs EHI mam: %.3f" % (len(H & M) / float(len(H | M))))
print("  jaccard newt vs EHI mam    : %.3f" % (len(N & M) / float(len(N | M))))
print()
print("  The last line is the key control: newts and mammals come from the SAME")
print("  study and the SAME pipeline, so any difference between them cannot be")
print("  a batch effect.")

key = sorted(H & N, key=lambda g: -herp_all[g])
print()
print("  GENERA IN BOTH AMPHIBIAN CATALOGS:")
print("    %-26s %8s %6s %8s %6s %6s" % ("genus", "herptile", "newt", "EHI-mam", "Yblut", "refs"))
for g in key:
    flag = "  <-- absent from both endotherm catalogs" if g not in endo else ""
    print("    %-26s %8d %6d %8d %6d %6d%s"
          % (g, herp_all[g], newt_all[g], ehi_mam.get(g, 0), yb.get(g, 0),
             ref.get(g, 0), flag))

both_amph_not_endo = [g for g in key if g not in endo]
mass_h = sum(herp_all[g] for g in both_amph_not_endo)
mass_n = sum(newt_all[g] for g in both_amph_not_endo)
print()
print("  IN BOTH AMPHIBIAN CATALOGS AND NEITHER ENDOTHERM CATALOG: %d genera"
      % len(both_amph_not_endo))
print("    holding %d herptile SGBs and %d newt genomes" % (mass_h, mass_n))
print("    This is the presence claim. It replaces an absence claim with a")
print("    positive one, replicated across two independent studies.")

print()
print("  herptile genera NOT in newts: %d" % len(H - N))
for g in sorted(H - N, key=lambda g: -herp_all[g])[:10]:
    print("    %-26s herptile %3d, refs %3d" % (g, herp_all[g], ref.get(g, 0)))

# ------------------------------------------- salamander vs frog
print()
print("=" * 74)
print("IS THE SHARING SALAMANDER-SPECIFIC?")
print("=" * 74)
print("  Newts are Caudata, like salamanders. If the shared genera are only")
print("  the salamander ones, the replication is within-order and does not")
print("  extend to frogs.")
for grp in ("Salamander", "Frog", "Toad", "Lizard"):
    s = set(herp_by_group.get(grp, {}))
    if not s:
        continue
    print("  %-12s %3d genera, %3d shared with newts (%.0f%%)"
          % (grp, len(s), len(s & N), 100.0 * len(s & N) / len(s)))

with open(OUT, "w") as f:
    f.write("genus\therptile_sgbs\tnewt_genomes\tehi_mammal\tyoungblut\tgtdb_refs\t"
            "in_both_amphibian\tin_any_endotherm\n")
    for g in sorted(set(herp_all) | set(newt_all) | endo | set(ref)):
        f.write("%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n"
                % (g, herp_all.get(g, 0), newt_all.get(g, 0), ehi_mam.get(g, 0),
                   yb.get(g, 0), ref.get(g, 0),
                   "yes" if (g in H and g in N) else "no",
                   "yes" if g in endo else "no"))
print()
print("wrote %s" % OUT)

# SENTINEL_END
