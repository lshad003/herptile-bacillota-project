#!/usr/bin/env python3
# Assessment of which SGBs are affected by excluding the laboratory-reared arm.
#
# Source: ruminococcaceae-agent/scripts/vivarium_sgb_impact.py
# Reads:  data/herptile_bacillota_A_HQ_manifest_with_source.tsv
#         data/sgb_manifest.tsv
# Writes: results/vivarium_sgb_impact.tsv
#
# Thirteen SGBs contain at least one genome from the laboratory-reared arm.
# Eleven consist only of such genomes and are lost entirely on exclusion. Two
# also contain genomes from wild-caught animals and are retained, but their
# representatives are drawn from the excluded arm and must be replaced.
#
# The distinction matters because taxonomy, completeness and contamination in
# the SGB manifest are taken from the representative rather than from the
# cluster, so a change of representative changes those fields.
import os, sys
from collections import defaultdict

BASE = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
MANI = os.path.join(BASE, "data/herptile_bacillota_A_HQ_manifest_with_source.tsv")
SGB  = os.path.join(BASE, "data/sgb_manifest.tsv")
OUT  = os.path.join(BASE, "results/vivarium_sgb_impact.tsv")

def die(m):
    sys.stderr.write("FATAL: " + m + "\n"); sys.exit(1)

if os.path.exists(OUT):
    die("output exists, refusing to overwrite: " + OUT)
for p in (MANI, SGB):
    if not os.path.isfile(p):
        die("missing " + p)

def read_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(head, l.rstrip("\n").split("\t")))
                for l in fh if l.strip()]
    return head, rows

mh, mrows = read_tsv(MANI)
print("manifest columns: " + ", ".join(mh))
if "source" not in mh:
    die("manifest lacks 'source'")
idcol = None
for c in ("genome", "mag", "bin", "genome_id", "mag_id", "bin_id"):
    if c in mh:
        idcol = c
        break
if idcol is None:
    die("cannot find a MAG id column in the manifest; columns are above")
print("using manifest MAG id column: " + idcol)
if len(mrows) != 2229:
    die("expected 2229 manifest rows, found %d" % len(mrows))

src_of = {}
for r in mrows:
    src_of[r[idcol]] = r["source"]
if len(src_of) != 2229:
    die("MAG ids are not unique: %d unique of 2229" % len(src_of))

sh, srows = read_tsv(SGB)
for c in ("sgb", "representative", "n_mags", "has_wild", "sources", "family"):
    if c not in sh:
        die("sgb_manifest lacks column " + c)
if len(srows) != 1171:
    die("expected 1171 SGB rows, found %d" % len(srows))

# Do the representatives join to manifest MAG ids at all?
hit = sum(1 for r in srows if r["representative"] in src_of)
print("")
print("representatives found in manifest by exact id: %d of 1171" % hit)
if hit < 1171:
    ex = [r["representative"] for r in srows if r["representative"] not in src_of][:5]
    print("  examples that did not join: " + ", ".join(ex))
    print("  example manifest ids       : " + ", ".join(list(src_of)[:5]))
    die("representative ids do not join to manifest MAG ids; stopping "
        "rather than guessing a transform")

# The `sources` column already records the per-SGB source set.
print("")
print("=== SGB SOURCE SETS, from sgb_manifest 'sources' ===")
by_src = defaultdict(int)
for r in srows:
    by_src[r["sources"]] += 1
for k in sorted(by_src, key=lambda x: -by_src[x]):
    print("  %-28s %d" % (k, by_src[k]))

viv = [r for r in srows if "VIVARIUM" in r["sources"]]
print("")
print("=== SGBs CONTAINING A VIVARIUM MAG: %d ===" % len(viv))
print("  %-14s %-8s %-7s %-9s %-22s %s"
      % ("sgb", "n_mags", "wild", "rep_src", "sources", "family"))
mixed = []
for r in sorted(viv, key=lambda x: x["sgb"]):
    rep_src = src_of.get(r["representative"], "NOT_FOUND")
    print("  %-14s %-8s %-7s %-9s %-22s %s"
          % (r["sgb"], r["n_mags"], r["has_wild"], rep_src,
             r["sources"], r["family"]))
    if r["has_wild"] == "yes":
        mixed.append((r, rep_src))

print("")
print("=== THE DECISION-RELEVANT CASES ===")
print("  vivarium-containing SGBs counted as WILD (has_wild=yes): %d" % len(mixed))
for r, rep_src in mixed:
    print("    %s  rep=%s  rep_source=%s  n_mags=%s  family=%s"
          % (r["sgb"], r["representative"], rep_src, r["n_mags"], r["family"]))
    if rep_src == "VIVARIUM":
        print("      ^ REPRESENTATIVE IS A VIVARIUM MAG. Removing the wood frog")
        print("        arm changes this SGB's representative genome, so any")
        print("        analysis keyed on it would need rebuilding.")

reps_viv = sum(1 for r in srows if src_of.get(r["representative"]) == "VIVARIUM")
print("")
print("  SGBs whose REPRESENTATIVE is a vivarium MAG, anywhere: %d" % reps_viv)

viv_only = [r for r in viv if r["has_wild"] != "yes"]
print("  vivarium-containing SGBs NOT counted as wild: %d" % len(viv_only))
print("")
print("IF THE COUNT OF WILD-COUNTED VIVARIUM SGBs IS 0, removing the wood frog")
print("arm cannot change the 718 wild SGBs and is a recount, not a rebuild.")

with open(OUT, "w") as out:
    out.write("sgb\trepresentative\trep_source\tn_mags\thas_wild\tsources\tfamily\n")
    for r in sorted(viv, key=lambda x: x["sgb"]):
        out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                  % (r["sgb"], r["representative"],
                     src_of.get(r["representative"], "NOT_FOUND"),
                     r["n_mags"], r["has_wild"], r["sources"], r["family"]))
print("")
print("WROTE: " + OUT)
# VIVARIUM_SGB_IMPACT_V1_20260808
