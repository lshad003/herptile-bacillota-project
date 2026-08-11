# Amphibian genomes are staged for biosynthetic gene cluster detection.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/stage_bgc_input.py
# Output: work/bgc/bgc_input_list.txt, bgc_manifest.tsv
# STAGE_BGC_INPUT_V1_20260805
import os
from collections import Counter, defaultdict

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
SGB = os.path.join(ROOT, "data/sgb_manifest.tsv")
HERP = os.path.join(ROOT, "data/herptile_bacillota_A_HQ_manifest_with_source.tsv")
HERP_DIR = os.path.join(ROOT, "results/drep_herptile_95ani_2229/dereplicated_genomes")
AMPH_MAN = os.path.join(ROOT, "results/ehi_amphibian_manifest.tsv")
AMPH_SUM = os.path.join(ROOT, "results/gtdbtk_ehi_amphibian_r220/gtdbtk.bac120.summary.tsv")

WORK = os.path.join(ROOT, "work/bgc")
LINKS = os.path.join(WORK, "genomes")
OUT_LIST = os.path.join(WORK, "bgc_input_list.txt")
OUT_MAN = os.path.join(WORK, "bgc_manifest.tsv")

AMPHIBIAN = {"Salamander", "Frog", "Toad"}
REPTILE = {"Lizard", "Turtle", "Tortoise", "Snake"}
EXT = (".fa", ".fna", ".fasta")


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


for d in (WORK, LINKS):
    if not os.path.isdir(d):
        os.makedirs(d)

taxon_map = {}
for r in read_tsv(HERP):
    t = r["host_taxon"].strip()
    if t:
        taxon_map[t] = r["animal_type"].strip()

# ---------------------------------------------------------- herptile wild
print("=" * 74)
print("ARM 1: WILD SGB REPRESENTATIVES (this study)")
print("=" * 74)
avail = {}
for e in os.scandir(HERP_DIR):
    if e.is_file() and e.name.endswith(EXT):
        stem = e.name
        for x in EXT:
            if stem.endswith(x):
                stem = stem[:-len(x)]
                break
        avail[stem] = e.path

entries = []
for r in read_tsv(SGB):
    if r["has_wild"].strip().lower() != "yes":
        continue
    groups = {taxon_map[t] for t in r["host_species"].split(";")
              if t.strip() and t.strip() in taxon_map}
    if groups & AMPHIBIAN and not groups & REPTILE:
        host = "amphibian"
    elif groups & REPTILE and not groups & AMPHIBIAN:
        host = "reptile"
    elif groups:
        host = "mixed"
    else:
        host = "unknown"
    rep = r["representative"].strip()
    p = avail.get(rep)
    entries.append(dict(arm="herptile", gid=rep, path=p or "",
                        family=r["family"].strip(), genus=r["genus"].strip(),
                        host=host, host_detail=r["host_species"].strip()))
found = sum(1 for e in entries if e["path"])
print("  wild SGBs: %d, FASTA located: %d, missing: %d"
      % (len(entries), found, len(entries) - found))
print("  by host: %s" % dict(Counter(e["host"] for e in entries)))
print("  by family (top 6): %s"
      % dict(Counter(e["family"] for e in entries).most_common(6)))

# ---------------------------------------------------------- ehi newts
print()
print("=" * 74)
print("ARM 2: EHI NEWT Bacillota_A")
print("=" * 74)
fam_of = {}
gen_of = {}
for r in read_tsv(AMPH_SUM):
    t = parse_tax(r["classification"])
    fam_of[r["user_genome"].strip()] = t.get("f", "")
    gen_of[r["user_genome"].strip()] = t.get("g", "")

n_amph = 0
for r in read_tsv(AMPH_MAN):
    g = r["genome_id"].strip()
    p = r["fasta"].strip()
    entries.append(dict(arm="ehi_amphibian", gid=g, path=p,
                        family=fam_of.get(g, ""), genus=gen_of.get(g, ""),
                        host="amphibian", host_detail=r["host_species"].strip()))
    n_amph += 1
print("  EHI newt genomes: %d" % n_amph)
print("  by host species: %s"
      % dict(Counter(e["host_detail"] for e in entries if e["arm"] == "ehi_amphibian")))
print("  by family (top 6): %s"
      % dict(Counter(e["family"] for e in entries
                     if e["arm"] == "ehi_amphibian").most_common(6)))
print()
print("  NOTE: Youngblut contributes no amphibian genomes. The 393 downloaded")
print("  explicitly excluded herptile and fish hosts, so it cannot appear in")
print("  an amphibian BGC survey.")

# ---------------------------------------------------------- resolve + link
print()
print("=" * 74)
print("RESOLVING AND LINKING")
print("=" * 74)
ok, missing = [], []
for e in entries:
    p = e["path"]
    if p and os.path.exists(p) and os.path.getsize(p) > 0:
        ok.append(e)
    else:
        missing.append(e)
print("  resolved: %d, missing: %d" % (len(ok), len(missing)))
for e in missing[:6]:
    print("     %-14s %-30s %s" % (e["arm"], e["gid"], e["path"] or "<no path>"))

# antiSMASH cannot read gzip, so gzipped inputs are flagged for the job to
# decompress into a scratch copy rather than being linked directly.
n_gz = sum(1 for e in ok if e["path"].endswith(".gz"))
print("  gzipped inputs: %d (the job decompresses these)" % n_gz)

n_link = 0
for e in ok:
    safe = "%s__%s" % (e["arm"], e["gid"].replace("/", "_"))
    e["safe"] = safe
    dst = os.path.join(LINKS, safe + (".fna.gz" if e["path"].endswith(".gz") else ".fna"))
    e["staged"] = dst
    if not os.path.exists(dst):
        os.symlink(e["path"], dst)
        n_link += 1
print("  symlinks created: %d in %s" % (n_link, LINKS))

with open(OUT_LIST, "w") as f:
    for e in ok:
        f.write("%s\n" % e["staged"])
with open(OUT_MAN, "w") as f:
    f.write("index\tsafe_name\tarm\tgenome\tfamily\tgenus\thost\thost_detail\tstaged_path\n")
    for i, e in enumerate(ok):
        f.write("%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                % (i, e["safe"], e["arm"], e["gid"], e["family"], e["genus"],
                   e["host"], e["host_detail"], e["staged"]))
print("  wrote %s (%d lines)" % (OUT_LIST, len(ok)))
print("  wrote %s" % OUT_MAN)

print()
print("=" * 74)
print("FINAL INPUT SET")
print("=" * 74)
print("  total genomes: %d" % len(ok))
print("  by arm : %s" % dict(Counter(e["arm"] for e in ok)))
print("  by host: %s" % dict(Counter(e["host"] for e in ok)))
amph_only = sum(1 for e in ok if e["host"] == "amphibian")
print("  amphibian-host genomes across both arms: %d" % amph_only)
print()
print("  SLURM array max is 1000, so %d needs the OFFSET pattern or two"
      % len(ok))
print("  submissions. Tasks 0-999 first, then an offset run for the rest.")

# SENTINEL_END
