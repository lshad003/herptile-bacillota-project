#!/usr/bin/env python3
# Rebuilding of the SGB manifest with the laboratory-reared arm excluded.
#
# Source: ruminococcaceae-agent/scripts/sgb_manifest_nowf.py
# Reads:  results/drep_herptile_95ani_2229/data_tables/Cdb.csv and Wdb.csv
#         data/herptile_bacillota_A_HQ_manifest_with_source.tsv
# Writes: data/sgb_manifest_nowf.tsv
#
# Exclusion is applied after dereplication rather than before, so cluster
# boundaries are unchanged and SGB identifiers remain comparable between the
# full and filtered manifests. Clusters left with no member are dropped;
# clusters whose representative is removed are assigned a new one from the
# remaining members.
#
# The unfiltered cluster count is reproduced before filtering is applied, as
# a check that the join between the dereplication tables and the MAG manifest
# is correct.
#
# strip_ext is imported from the original manifest builder rather than
# reimplemented, so genome identifiers are transformed identically in both.
import os, sys, importlib.util
from collections import defaultdict

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
DREP = B + "/results/drep_herptile_95ani_2229/data_tables"
CDB = DREP + "/Cdb.csv"
WDB = DREP + "/Wdb.csv"
MAN = B + "/data/herptile_bacillota_A_HQ_manifest_with_source.tsv"
ORIG = B + "/data/sgb_manifest.tsv"
OUT = B + "/data/sgb_manifest_nowf.tsv"
SRC_BUILDER = B + "/scripts/build_sgb_manifest.py"

DROP_SOURCE = "VIVARIUM"

def die(m):
    sys.stderr.write("FATAL: " + m + "\n"); sys.exit(1)

if os.path.exists(OUT):
    die("output exists, refusing to overwrite: " + OUT)
for p in (CDB, MAN, ORIG, SRC_BUILDER):
    if not os.path.isfile(p):
        die("missing " + p)

# Reuse strip_ext from the original builder rather than reimplementing it,
# so cluster-to-genome joins are identical between the two manifests.
src = open(SRC_BUILDER).read()
ns = {}
for line in src.split("\n"):
    if line.startswith("def strip_ext"):
        break
else:
    die("strip_ext not found in " + SRC_BUILDER)
start = src.index("def strip_ext")
tail = src[start:]
end = len(tail)
for i, line in enumerate(tail.split("\n")[1:], start=1):
    if line and not line[0].isspace():
        end = sum(len(x) + 1 for x in tail.split("\n")[:i])
        break
exec(compile(tail[:end], "strip_ext", "exec"), ns)
strip_ext = ns.get("strip_ext")
if strip_ext is None:
    die("could not extract strip_ext")
print("reusing strip_ext from build_sgb_manifest.py")

def read_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return head, [dict(zip(head, l.rstrip("\n").split("\t")))
                      for l in fh if l.strip()]

mh, mrows = read_tsv(MAN)
if len(mrows) != 2229:
    die("expected 2229 manifest rows, found %d" % len(mrows))
mag = {}
for r in mrows:
    try:
        comp = float(r["completeness"]); cont = float(r["contamination"])
    except ValueError:
        die("unparseable completeness/contamination for " + r["bin_id"])
    mag[r["bin_id"]] = dict(src=r["source"], comp=comp, cont=cont,
                            fam=r.get("family", ""), tax=r.get("taxonomy", ""),
                            base=r.get("sample_id_base", ""))
if len(mag) != 2229:
    die("bin_id not unique")

drop = set(g for g, d in mag.items() if d["src"] == DROP_SOURCE)
print("MAGs dropped (%s): %d of 2229" % (DROP_SOURCE, len(drop)))
print("MAGs retained: %d" % (2229 - len(drop)))

clu = {}
with open(CDB) as fh:
    ch = fh.readline().rstrip("\n").replace('"', "").split(",")
    cg, cs = ch.index("genome"), ch.index("secondary_cluster")
    for line in fh:
        p = line.rstrip("\n").replace('"', "").split(",")
        if len(p) > max(cg, cs):
            clu[strip_ext(p[cg])] = p[cs]
print("genomes in Cdb: %d | clusters: %d" % (len(clu), len(set(clu.values()))))

winner = {}
if os.path.exists(WDB):
    with open(WDB) as fh:
        wh = fh.readline().rstrip("\n").replace('"', "").split(",")
        wg = wh.index("genome")
        wc = wh.index("cluster") if "cluster" in wh else None
        for line in fh:
            p = line.rstrip("\n").replace('"', "").split(",")
            if len(p) <= wg:
                continue
            g = strip_ext(p[wg])
            c = p[wc] if wc is not None and len(p) > wc else clu.get(g)
            if c:
                winner[c] = g
    print("Wdb winners: %d" % len(winner))
else:
    die("no Wdb.csv; refusing to re-pick every representative by fallback")

full = defaultdict(list)
for g, c in clu.items():
    if g in mag:
        full[c].append(g)
if len(full) != 1171:
    die("unfiltered rebuild gives %d clusters, expected 1171" % len(full))
print("unfiltered clusters reproduce: 1171")

kept = defaultdict(list)
for c, mem in full.items():
    m = [g for g in mem if g not in drop]
    if m:
        kept[c] = m

lost = sorted(set(full) - set(kept))
print("")
print("=== SGBs LOST ENTIRELY: %d ===" % len(lost))
for c in lost:
    mem = full[c]
    fam = mag[winner.get(c, mem[0])]["fam"] if winner.get(c, mem[0]) in mag else "?"
    print("  %-10s n_mags=%-3d family=%s" % (c, len(mem), fam))

rows = []
repicked = []
taxchange = []
for c, mem in sorted(kept.items()):
    rep = winner.get(c)
    if rep not in mag or rep in drop:
        old = rep
        rep = max(mem, key=lambda g: mag[g]["comp"] - 5.0 * mag[g]["cont"])
        repicked.append((c, old, rep, len(full[c]), len(mem)))
        if old in mag and mag[old]["fam"] != mag[rep]["fam"]:
            taxchange.append((c, old, mag[old]["fam"], rep, mag[rep]["fam"]))
    r = mag[rep]
    srcs = sorted(set(mag[g]["src"] for g in mem))
    rows.append(dict(sgb=c, rep=rep, n=len(mem), fam=r["fam"],
                     comp=r["comp"], cont=r["cont"],
                     src=",".join(srcs), wild=("WILD" in srcs)))

print("")
print("=== REPRESENTATIVES RE-PICKED: %d ===" % len(repicked))
for c, old, new, n_before, n_after in repicked:
    print("  %-10s %s -> %s  (n %d -> %d)" % (c, old, new, n_before, n_after))
    print("             old fam=%s  new fam=%s"
          % (mag[old]["fam"] if old in mag else "?", mag[new]["fam"]))

if taxchange:
    print("")
    print("  WARNING: FAMILY ASSIGNMENT CHANGED IN %d SGB(s)." % len(taxchange))
    for c, old, of, new, nf in taxchange:
        print("    %s: %s (%s) -> %s (%s)" % (c, old, of, new, nf))
    print("  Re-run GTDB-Tk on the new representative before using these rows.")

oh, orows = read_tsv(ORIG)
o_wild = sum(1 for r in orows if r["has_wild"] == "yes")
o_rum = sum(1 for r in orows if r["family"] == "Ruminococcaceae")
n_wild = sum(1 for r in rows if r["wild"])
n_rum = sum(1 for r in rows if r["fam"] == "Ruminococcaceae")

print("")
print("=== COUNTS, WITH AND WITHOUT THE WOOD FROG ARM ===")
print("  %-28s %8s %8s %8s" % ("", "with", "without", "delta"))
print("  %-28s %8d %8d %8d" % ("MAGs", 2229, 2229 - len(drop), -len(drop)))
print("  %-28s %8d %8d %8d" % ("SGBs", len(orows), len(rows), len(rows) - len(orows)))
print("  %-28s %8d %8d %8d" % ("wild-containing SGBs", o_wild, n_wild, n_wild - o_wild))
print("  %-28s %8d %8d %8d" % ("Ruminococcaceae SGBs", o_rum, n_rum, n_rum - o_rum))
print("")
print("  NOTE family here comes from the MAG manifest taxonomy, not from a")
print("  fresh GTDB-Tk run. Any re-picked representative needs GTDB-Tk,")
print("  CheckM, GUNC and rRNA/tRNA re-run before its row is quotable.")

with open(OUT, "w") as f:
    f.write("sgb\trepresentative\tn_mags\tfamily\trep_completeness\t"
            "rep_contamination\tsources\thas_wild\trepicked\n")
    rp = set(c for c, _, _, _, _ in repicked)
    for r in rows:
        f.write("%s\t%s\t%d\t%s\t%.2f\t%.2f\t%s\t%s\t%s\n"
                % (r["sgb"], r["rep"], r["n"], r["fam"], r["comp"], r["cont"],
                   r["src"], "yes" if r["wild"] else "no",
                   "yes" if r["sgb"] in rp else "no"))
print("")
print("WROTE: " + OUT)
print("This file is DELIBERATELY REDUCED: it carries the columns needed to")
print("recount, not the full 18-column schema. It is for deciding, not for")
print("replacing data/sgb_manifest.tsv.")
# SGB_MANIFEST_NOWF_V1_20260808
