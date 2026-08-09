#!/usr/bin/env python3
# Builds the three supplementary tables for the catalog step.
# Source: ruminococcaceae-agent/scripts/make_step1_tables_v2.py
# Output: tables/TableS1_strict_chimerism.tsv
#         tables/TableS2_catalog_by_family.tsv
#         tables/TableS3_feature_recovery.tsv
#
# Feature recovery denominators differ between rRNA and tRNA where a tool
# failed on a genome: barrnap succeeded on all, tRNAscan-SE failed on one
# reference.
import os, sys
from collections import Counter

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
MAN = B + "/data/herptile_bacillota_A_HQ_manifest_with_source.tsv"
SGB = B + "/data/sgb_manifest.tsv"
GUNC = B + "/results/gunc_audit_by_sgb.tsv"
CAT = B + "/results/rrna_trna_catalog_per_genome_merged.tsv"
REF = B + "/results/rrna_trna_refs_per_genome.tsv"

S1 = B + "/results/TableS1_strict_chimerism_v2.tsv"
S2 = B + "/results/TableS2_catalog_by_family_v2.tsv"
S3 = B + "/results/TableS3_feature_recovery_v3.tsv"

def die(m):
    sys.stderr.write("FATAL: " + m + "\n"); sys.exit(1)

for p in (S1, S2, S3):
    if os.path.exists(p):
        die("output exists, refusing to overwrite: " + p)
for p in (MAN, SGB, GUNC, CAT, REF):
    if not os.path.isfile(p):
        die("missing " + p)

def read_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return head, [dict(zip(head, l.rstrip("\n").split("\t")))
                      for l in fh if l.strip()]

mh, mrows = read_tsv(MAN)
sh, srows = read_tsv(SGB)
gh, grows = read_tsv(GUNC)
ch, crows = read_tsv(CAT)
rh, rrows = read_tsv(REF)

if len(mrows) != 2229: die("expected 2229 MAG rows, found %d" % len(mrows))
if len(srows) != 1171: die("expected 1171 SGB rows, found %d" % len(srows))
if len(grows) != 1171: die("expected 1171 GUNC rows, found %d" % len(grows))

# -------------------------------------------------- S1 strict chimerism
strict = [r for r in grows if r["strict_chimera"] == "yes"]
print("=== S1: strict chimerism calls ===")
print("  rows: %d" % len(strict))
for k, v in Counter(r["family"] for r in strict).most_common():
    print("    %-24s %d" % (k, v))
print("  passing GUNC's own criterion: %d"
      % sum(1 for r in strict if r["pass_gunc"] == "pass"))

cols1 = ["genome", "sgb", "family", "genus", "taxonomic_level",
         "clade_separation_score", "reference_representation_score",
         "checkm_completeness", "checkm_contamination", "n_contigs",
         "pass_gunc"]
with open(S1, "w") as f:
    f.write("\t".join(cols1) + "\n")
    for r in sorted(strict, key=lambda x: -float(x["clade_separation_score"])):
        f.write("\t".join(r.get(c, "") for c in cols1) + "\n")

# -------------------------------------------------- S2 composition by family
mag_f = Counter(r["family"] for r in mrows)
sgb_f = Counter(r["family"] for r in srows)
wild_f = Counter(r["family"] for r in srows if r["has_wild"] == "yes")
sing_f = Counter(r["family"] for r in srows if r["n_mags"] == "1")

rows2 = [((fm if fm else "(unassigned)"), mag_f.get(fm, 0), sgb_f[fm],
          wild_f.get(fm, 0), sing_f.get(fm, 0))
         for fm in sorted(sgb_f, key=lambda f: -mag_f.get(f, 0))]

print("")
print("=== S2: catalog composition by family ===")
print("  %-26s %7s %7s %7s %7s" % ("family", "MAGs", "SGBs", "wild", "single"))
for r in rows2[:10]:
    print("  %-26s %7d %7d %7d %7d" % r)
print("  ... %d families in total" % len(rows2))
print("  %-26s %7d %7d %7d %7d"
      % ("TOTAL", sum(x[1] for x in rows2), sum(x[2] for x in rows2),
         sum(x[3] for x in rows2), sum(x[4] for x in rows2)))

with open(S2, "w") as f:
    f.write("family\tn_mags\tn_sgbs\tn_wild_sgbs\tn_singleton_sgbs\n")
    for r in rows2:
        f.write("%s\t%d\t%d\t%d\t%d\n" % r)

# -------------------------------------------------- S3 feature recovery
# Denominators differ by feature where a tool failed on a genome: barrnap
# succeeded on every genome in both arms, tRNAscan-SE failed on one reference.
rum = set(r["representative"] for r in srows if r["family"] == "Ruminococcaceae")

def feat(rows, label, subset=None):
    rs = [r for r in rows if subset is None or r["genome"] in subset]
    ok_b = [r for r in rs if r.get("barrnap_ok", "1") == "1"]
    ok_t = [r for r in rs if r.get("trnascan_ok", "1") == "1"]
    ok_x = [r for r in rs if r.get("barrnap_ok", "1") == "1"
            and r.get("trnascan_ok", "1") == "1"]
    if not ok_b or not ok_t:
        die("no parseable rows for " + label)
    p = lambda sub, fn: 100.0 * sum(1 for r in sub if fn(r)) / len(sub)
    return dict(arm=label, n_rrna=len(ok_b), n_trna=len(ok_t),
        pct_5S=p(ok_b, lambda r: r["s5"] == "complete"),
        pct_16S=p(ok_b, lambda r: r["s16"] == "complete"),
        pct_23S=p(ok_b, lambda r: r["s23"] == "complete"),
        pct_all_three=p(ok_b, lambda r: r["all_three_complete"] == "1"),
        pct_trna=p(ok_t, lambda r: r["aa18"] == "1"),
        pct_both=p(ok_x, lambda r: r["all_three_complete"] == "1"
                   and r["aa18"] == "1"))

rows3 = [feat(crows, "Catalog, all SGB representatives"),
         feat(crows, "Catalog, Ruminococcaceae", rum),
         feat(rrows, "GTDB r220 Ruminococcaceae references")]

print("")
print("=== S3: rRNA and tRNA recovery ===")
print("  %-38s %6s %6s %7s %7s %7s %8s %7s"
      % ("arm", "n rRNA", "n tRNA", "5S", "16S", "23S", "all 3", ">=18AA"))
for d in rows3:
    print("  %-38s %6d %6d %6.1f%% %6.1f%% %6.1f%% %7.1f%% %6.1f%%"
          % (d["arm"], d["n_rrna"], d["n_trna"], d["pct_5S"], d["pct_16S"],
             d["pct_23S"], d["pct_all_three"], d["pct_trna"]))

cols3 = ["arm", "n_rrna", "n_trna", "pct_5S", "pct_16S", "pct_23S",
         "pct_all_three", "pct_trna", "pct_both"]
with open(S3, "w") as f:
    f.write("\t".join(cols3) + "\n")
    for d in rows3:
        f.write("\t".join(("%s" % d[c]) if c == "arm" else
                          ("%d" % d[c]) if c.startswith("n_") else
                          ("%.1f" % d[c]) for c in cols3) + "\n")

print("")
print("WROTE:")
for p in (S1, S2, S3):
    print("  " + p)
# MAKE_STEP1_TABLES_V2_20260809
