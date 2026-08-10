#!/usr/bin/env python3
# Supplementary table S7 is built for section 3.4.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/make_step4_table.py
# Output: tables/TableS7_read_recovery_testability.tsv
"""
Supplementary table S7 for section 3.4.
Wild catalog genera with their testability class, plus captive-only genera
appearing in the recovery join. Refuses to overwrite an existing file.
"""
import os
import sys

AG = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
RP = "/bigdata/stajichlab/lshad003/herptile-bacillota-project"

RECIP = os.path.join(AG, "results/gtdb_ncbi_reciprocal.tsv")
JOIN = os.path.join(AG, "results/bracken_reciprocal_recovery_v2.tsv")
OUT = os.path.join(RP, "tables/TableS7_read_recovery_testability.tsv")

if os.path.exists(OUT):
    sys.exit("REFUSING TO OVERWRITE %s" % OUT)

recip = {}
with open(RECIP) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ix = {n: i for i, n in enumerate(hdr)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        recip[f[ix["gtdb_genus"]]] = {
            "family": f[ix["gtdb_family"]],
            "wild": f[ix["in_wild_catalog"]].strip().lower(),
            "ncbi": f[ix["best_ncbi_genus"]],
            "fwd": f[ix["forward_share"]],
            "rev": f[ix["reverse_share"]],
            "verdict": f[ix["verdict_t70"]],
        }

joined = {}
with open(JOIN) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    jx = {n: i for i, n in enumerate(hdr)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        joined[f[jx["gtdb_genus"]]] = {
            "bucket": f[jx["bucket"]],
            "n_samples": f[jx["n_wild_samples_detected"]],
            "mean": f[jx["mean_fraction_wild"]],
            "max": f[jx["max_fraction_wild"]],
            "wild_sgbs": f[jx["wild_sgbs"]],
            "captive_sgbs": f[jx["captive_sgbs"]],
        }

wildset = set(g for g, d in recip.items() if d["wild"] == "yes")
keep = sorted(wildset | set(joined))
print("wild catalog genera: %d" % len(wildset))
print("genera in the recovery join: %d" % len(joined))
print("captive-only genera added from the join: %d" % len(set(joined) - wildset))
print("table rows: %d" % len(keep))

missing_recip = [g for g in keep if g not in recip]
if missing_recip:
    print("NOT IN THE RECIPROCAL FILE AT ALL: %d" % len(missing_recip))
    for m in missing_recip[:10]:
        print("  %s" % m)

print("")
print("TESTABILITY CLASSES, WILD CATALOG GENERA ONLY")
cls = {}
for g in wildset:
    v = recip[g]["verdict"]
    cls[v] = cls.get(v, 0) + 1
for k in sorted(cls, key=lambda x: -cls[x]):
    print("  %-34s %d" % (k, cls[k]))

print("")
print("BUCKETS, SPLIT BY WHETHER THE GENUS IS IN THE WILD CATALOG")
bk = {}
for g, j in joined.items():
    key = (j["bucket"], "wild" if g in wildset else "captive_only")
    bk[key] = bk.get(key, 0) + 1
for k in sorted(bk, key=lambda x: -bk[x]):
    print("  %-32s %-14s %d" % (k[0], k[1], bk[k]))

rows = []
for g in keep:
    d = recip.get(g, {"family": "", "wild": "", "ncbi": "",
                      "fwd": "", "rev": "", "verdict": ""})
    j = joined.get(g)
    rows.append([
        g, d["family"], "yes" if g in wildset else "no", d["ncbi"],
        d["fwd"], d["rev"], d["verdict"],
        j["bucket"] if j else "not_testable",
        j["n_samples"] if j else "",
        j["mean"] if j else "",
        j["max"] if j else "",
        j["wild_sgbs"] if j else "",
        j["captive_sgbs"] if j else "",
    ])

with open(OUT, "w") as fh:
    fh.write("gtdb_genus\tgtdb_family\tin_wild_catalog\tbest_ncbi_genus\t"
             "forward_share\treverse_share\ttestability_t70\tbucket\t"
             "n_wild_samples_detected\tmean_fraction_wild\tmax_fraction_wild\t"
             "wild_sgbs\tcaptive_sgbs\n")
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")

print("")
print("WROTE %s, %d rows" % (OUT, len(rows)))
print("")
print("CAPTION: read fractions are given only for genera passing the")
print("reciprocal 70 percent criterion in both directions. Genera with")
print("in_wild_catalog = no are recovered only from captive animals and are")
print("included because they carry the read-versus-recovery contrast. A")
print("bucket of not_testable means the genus could not be tested, which is")
print("not the same as not being detected. Read classifications name a taxon")
print("rather than demonstrating its presence.")
# MAKE_STEP4_TABLE_V2
