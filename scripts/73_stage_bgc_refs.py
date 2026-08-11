# GTDB reference genomes are staged for biosynthetic gene cluster detection.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/stage_bgc_refs.py
# Output: work/bgc_refs/bgc_ref_input_list.txt, bgc_ref_manifest.tsv
# STAGE_BGC_REFS_V1_20260805
import os
from collections import Counter

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
REFS = os.path.join(ROOT, "data/tasks/gtdb_ruminococcaceae_paths.tsv")
GTDB_TAX = "/srv/projects/db/gtdbtk/220/taxonomy/gtdb_taxonomy.tsv"
WORK = os.path.join(ROOT, "work/bgc_refs")
LINKS = os.path.join(WORK, "genomes")
OUT_LIST = os.path.join(WORK, "bgc_ref_input_list.txt")
OUT_MAN = os.path.join(WORK, "bgc_ref_manifest.tsv")


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

genus = {}
with open(GTDB_TAX) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        t = parse_tax(f[1])
        if t.get("f", "") == "Ruminococcaceae":
            genus[f[0].strip().replace("GB_", "").replace("RS_", "")] = t.get("g", "")

refs = read_tsv(REFS)
print("=" * 74)
print("GTDB r220 RUMINOCOCCACEAE REFERENCES")
print("=" * 74)
print("  rows in paths file: %d" % len(refs))
print()
print("  The amphibian antiSMASH run has no comparison arm: 1,155 genomes")
print("  compared to nothing. This adds the reference genomes so BGC content")
print("  can be contrasted. Note that references are 92%% MAGs themselves")
print("  (CHATINDEX R1), so this is a MAG-to-MAG comparison and fragmentation")
print("  affects both arms, which is what makes it fair.")

ok, missing = [], []
for r in refs:
    p = r["path"].strip()
    if p and os.path.exists(p) and os.path.getsize(p) > 0:
        ok.append(r)
    else:
        missing.append(r)
print()
print("  resolved: %d, missing: %d" % (len(ok), len(missing)))
for r in missing[:5]:
    print("     %s %s" % (r["accession"], r["path"]))

n_link = 0
entries = []
for r in ok:
    acc = r["accession"].strip()
    safe = "gtdbref__" + acc
    src = r["path"].strip()
    dst = os.path.join(LINKS, safe + (".fna.gz" if src.endswith(".gz") else ".fna"))
    if not os.path.exists(dst):
        os.symlink(src, dst)
        n_link += 1
    entries.append(dict(safe=safe, acc=acc, genus=genus.get(acc, ""),
                        staged=dst, src=src,
                        in_filtered=r.get("in_host_filtered_set", "")))
print("  symlinks created: %d" % n_link)
print("  gzipped inputs: %d (the job decompresses these)"
      % sum(1 for e in entries if e["staged"].endswith(".gz")))

with open(OUT_LIST, "w") as f:
    for e in entries:
        f.write("%s\n" % e["staged"])
with open(OUT_MAN, "w") as f:
    f.write("index\tsafe_name\taccession\tgenus\tstaged_path\tsource_path\n")
    for i, e in enumerate(entries):
        f.write("%d\t%s\t%s\t%s\t%s\t%s\n"
                % (i, e["safe"], e["acc"], e["genus"], e["staged"], e["src"]))
print()
print("  wrote %s (%d lines)" % (OUT_LIST, len(entries)))
print("  wrote %s" % OUT_MAN)

print()
print("  top genera in the reference set:")
for g, n in Counter(e["genus"] for e in entries).most_common(10):
    print("     %-26s %4d" % (g if g else "<none>", n))
print()
print("  SLURM array max is 1000 and this is %d genomes, so the job is"
      % len(entries))
print("  submitted twice with OFFSET, as for the amphibian run.")

# SENTINEL_END
