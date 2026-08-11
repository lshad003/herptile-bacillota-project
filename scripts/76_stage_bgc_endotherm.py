# Endotherm comparison genomes are staged, selected on GTDB-Tk classification.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/stage_bgc_endotherm.py
# Output: work/bgc_endo/ input list and manifest
# STAGE_BGC_ENDOTHERM_V1_20260806
# Builds the input list and manifest for jobs/52_antismash_endotherm.sh.
# 280 EHI mammal/bird Ruminococcaceae + 48 Youngblut Ruminococcaceae,
# selected on column 2 (classification) of each GTDB-Tk summary.
#
# Genomes are NOT copied. The list points at the originals:
#   EHI       ch3-chitin-evolution/data/ehi_2025/mags/nonherptile_fa/<id>.fa.gz
#   Youngblut data/youngblut/genomes/<id>.contigs.fa.gz
# NOTE the Youngblut filenames carry ".contigs" before ".fa.gz"; the batchfile
# is the authoritative path-to-id mapping and is used rather than guessed.
#
# Safe names are prefixed by arm so the parser can split them later, matching
# the convention in work/bgc (herptile__, ehi_amphibian__) and work/bgc_refs
# (gtdbref__).

import os, sys

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
EHI_SUM = ROOT + "/results/gtdbtk_ehi_r220_classify/gtdbtk.bac120.summary.tsv"
YB_SUM = ROOT + "/results/gtdbtk_youngblut_r220/gtdbtk.bac120.summary.tsv"
EHI_MAN = ROOT + "/results/ehi_nonherptile_manifest.tsv"
YB_BATCH = ROOT + "/data/youngblut/youngblut_gtdbtk_batchfile.tsv"
EHI_DIR = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution/data/ehi_2025/mags/nonherptile_fa"
WORK = ROOT + "/work/bgc_endo"
LIST = WORK + "/bgc_endo_input_list.txt"
MAN = WORK + "/bgc_endo_manifest.tsv"

for p in (LIST, MAN):
    if os.path.exists(p):
        raise SystemExit("REFUSING TO OVERWRITE %s, move it first" % p)
for p in (EHI_SUM, YB_SUM, EHI_MAN, YB_BATCH):
    if not os.path.exists(p):
        sys.exit("MISSING: %s" % p)
if not os.path.isdir(EHI_DIR):
    sys.exit("MISSING: %s" % EHI_DIR)
os.makedirs(WORK, exist_ok=True)


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


ehi_meta = {r["genome_id"]: r for r in read_tsv(EHI_MAN)}
print("EHI manifest rows: %d" % len(ehi_meta))

yb_path = {}
with open(YB_BATCH) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            yb_path[f[1].strip()] = f[0].strip()
print("Youngblut batchfile entries: %d" % len(yb_path))

entries = []
missing = []

for r in read_tsv(EHI_SUM):
    t = parse_tax(r["classification"])
    if t.get("f", "") != "Ruminococcaceae":
        continue
    gid = r["user_genome"].strip()
    fa = os.path.join(EHI_DIR, gid + ".fa.gz")
    if not os.path.exists(fa):
        missing.append(("ehi_mammal", gid, fa))
        continue
    m = ehi_meta.get(gid, {})
    entries.append(dict(arm="ehi_mammal", safe="ehimam__" + gid, gid=gid,
                        path=fa, genus=t.get("g", ""),
                        host_species=m.get("host_species", ""),
                        host_class=m.get("host_class", ""),
                        host_order=m.get("host_order", ""),
                        completeness=m.get("completeness", ""),
                        contigs=m.get("contigs", "")))

for r in read_tsv(YB_SUM):
    t = parse_tax(r["classification"])
    if t.get("f", "") != "Ruminococcaceae":
        continue
    gid = r["user_genome"].strip()
    fa = yb_path.get(gid)
    if fa is None or not os.path.exists(fa):
        missing.append(("youngblut", gid, fa or "no batchfile entry"))
        continue
    entries.append(dict(arm="youngblut", safe="youngblut__" + gid, gid=gid,
                        path=fa, genus=t.get("g", ""),
                        host_species="", host_class="", host_order="",
                        completeness="", contigs=""))

n_ehi = sum(1 for e in entries if e["arm"] == "ehi_mammal")
n_yb = sum(1 for e in entries if e["arm"] == "youngblut")
print()
print("staged: %d total" % len(entries))
print("  ehi_mammal : %d  (expected 280)" % n_ehi)
print("  youngblut  : %d  (expected 48)" % n_yb)
print("missing: %d" % len(missing))
for a, g, p in missing[:10]:
    print("   %s %s %s" % (a, g, p))

if n_ehi != 280:
    sys.exit("EXPECTED 280 EHI Ruminococcaceae, got %d. Stopping." % n_ehi)
if n_yb != 48:
    sys.exit("EXPECTED 48 Youngblut Ruminococcaceae, got %d. Stopping." % n_yb)

hosts = {}
for e in entries:
    if e["arm"] == "ehi_mammal" and e["host_order"]:
        hosts[e["host_order"]] = hosts.get(e["host_order"], 0) + 1
print()
print("EHI mammal host orders:")
for k in sorted(hosts, key=lambda x: -hosts[x]):
    print("  %-24s %4d" % (k, hosts[k]))

entries.sort(key=lambda e: e["safe"])
with open(LIST, "w") as f:
    for e in entries:
        f.write("%s\t%s\n" % (e["path"], e["safe"]))
with open(MAN, "w") as f:
    f.write("index\tsafe_name\tarm\tgenome_id\tgenus\thost_species\t"
            "host_class\thost_order\tcompleteness\tcontigs\tstaged_path\n")
    for i, e in enumerate(entries):
        f.write("%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                % (i, e["safe"], e["arm"], e["gid"], e["genus"],
                   e["host_species"], e["host_class"], e["host_order"],
                   e["completeness"], e["contigs"], e["path"]))
print()
print("wrote %s  (%d lines)" % (LIST, len(entries)))
print("wrote %s" % MAN)
print()
print("ARRAY SIZE: %d genomes needs --array=0-%d" % (len(entries), len(entries) - 1))
print("jobs/52_antismash_endotherm.sh is written for 0-327. If the count")
print("above is not 328, FIX THE ARRAY BEFORE SUBMITTING.")
print("STAGE_BGC_ENDOTHERM_V1_20260806_COMPLETE")
