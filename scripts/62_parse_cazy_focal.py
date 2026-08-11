# Domain hits are filtered at the dbCAN published cutoffs and reduced to a family matrix.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/parse_cazy_focal.py
# Output: results/cazy_focal_hits.tsv, results/cazy_focal_family_matrix.tsv
import os, sys, csv, re
from collections import Counter, defaultdict

BASE = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
WORK = os.path.join(BASE, "work/cazy_focal")
DOMD = os.path.join(WORK, "domtbl")
LIST = os.path.join(WORK, "cazy_targets.tsv")
OUTH = os.path.join(BASE, "results/cazy_focal_hits.tsv")
OUTM = os.path.join(BASE, "results/cazy_focal_family_matrix.tsv")

EVAL = 1e-15
COV  = 0.35

def die(m):
    print("")
    print("!" * 72)
    print("FAILED: " + m)
    print("!" * 72)
    sys.exit(1)

def fam(name):
    n = name
    if n.endswith(".hmm"):
        n = n[:-4]
    return n

def cls(f):
    m = re.match(r"^(GH|GT|PL|CE|AA|CBM)", f)
    return m.group(1) if m else "other"

print("=" * 72)
print("STEP 1  TARGETS")
print("=" * 72)
if not os.path.exists(LIST):
    die("missing " + LIST)
rows = list(csv.DictReader(open(LIST), delimiter="\t"))
meta = {r["genome_id"]: r for r in rows}
print("  targets: %d" % len(rows))
print("  by set  : %s" % dict(Counter(r["set"] for r in rows)))
print("  by genus: %s" % dict(Counter(r["genus"] for r in rows).most_common(6)))
print("  by arm  : %s" % dict(Counter(r["arm"] for r in rows)))

files = {e.name[:-7]: e.path for e in os.scandir(DOMD) if e.name.endswith(".domtbl")}
print("  domtbl files: %d" % len(files))
missing = [g for g in meta if g not in files]
if missing:
    print("  MISSING: %s" % missing[:8])
    die("not every target has a domtbl")

print("")
print("=" * 72)
print("STEP 2  PARSE AND FILTER")
print("=" * 72)
print("  dbCAN thresholds: domain i-Evalue < %g, HMM coverage > %.2f" % (EVAL, COV))
raw = kept = fail_e = fail_c = 0
hits = []
percount = {}
for g, p in files.items():
    n = 0
    with open(p) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 23:
                continue
            raw += 1
            prot, hmm = f[0], fam(f[3])
            hlen = int(f[5])
            ie = float(f[12])
            hfrom, hto = int(f[15]), int(f[16])
            cov = (hto - hfrom + 1) / float(hlen) if hlen else 0.0
            if ie >= EVAL:
                fail_e += 1
                continue
            if cov <= COV:
                fail_c += 1
                continue
            kept += 1
            n += 1
            hits.append((g, prot, hmm, cls(hmm), ie, cov))
    percount[g] = n
print("  raw domain rows : %d" % raw)
print("  removed by E    : %d" % fail_e)
print("  removed by cov  : %d" % fail_c)
print("  kept            : %d" % kept)
if kept == 0:
    die("no hits survived, thresholds or parsing are wrong")

byfam = Counter(h[2] for h in hits)
bycls = Counter(h[3] for h in hits)
print("  distinct families: %d" % len(byfam))
print("  by class: %s" % dict(bycls.most_common()))
print("  top families: %s" % dict(byfam.most_common(10)))

print("")
print("=" * 72)
print("STEP 3  PER GENOME")
print("=" * 72)
gm = defaultdict(lambda: Counter())
for g, prot, hmm, c, ie, cov in hits:
    gm[g][hmm] += 1
counts = {g: sum(gm[g].values()) for g in meta}
byset = defaultdict(list)
for g in meta:
    byset[meta[g]["set"]].append(counts.get(g, 0))
for s in sorted(byset):
    v = sorted(byset[s])
    print("  %-18s n=%3d  median %4d  min %3d  max %4d"
          % (s, len(v), v[len(v)//2], v[0], v[-1]))

print("")
print("  FOCAL 125, BY ARM (this is the tested contrast):")
byarm = defaultdict(list)
for g in meta:
    if meta[g]["set"] == "focal125":
        byarm[meta[g]["arm"]].append(counts.get(g, 0))
for a in sorted(byarm):
    v = sorted(byarm[a])
    print("    %-12s n=%3d  median %4d  min %3d  max %4d"
          % (a, len(v), v[len(v)//2], v[0], v[-1]))
print("    NOT A TEST. Genome completeness differs by arm (amphibian 87.6 vs")
print("    reference 94.0 CheckM2), so a raw count difference is expected")
print("    from completeness alone. Per-family prevalence with a completeness")
print("    model is the test, not this table.")

print("")
print("  UNASSIGNED CLADES, descriptive only, no matched comparison group:")
byclade = defaultdict(list)
for g in meta:
    if meta[g]["set"] == "unassigned_clade":
        byclade[meta[g]["clade"]].append(counts.get(g, 0))
for c in sorted(byclade, key=lambda x: (x == "singleton", x)):
    v = sorted(byclade[c])
    print("    clade %-10s n=%2d  median %4d  min %3d  max %4d"
          % (c, len(v), v[len(v)//2], v[0], v[-1]))

if os.path.exists(OUTH):
    print("")
    print("  NOT overwriting existing " + OUTH)
else:
    with open(OUTH, "w") as fh:
        fh.write("genome_id\tset\tgenus\tarm\tclade\tprotein\tfamily\tclass\tievalue\thmm_coverage\n")
        for g, prot, hmm, c, ie, cov in hits:
            r = meta[g]
            fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.3g\t%.3f\n"
                     % (g, r["set"], r["genus"], r["arm"], r["clade"], prot, hmm, c, ie, cov))
    print("")
    print("  wrote %s (%d rows)" % (OUTH, kept))

fams = sorted(byfam)
if os.path.exists(OUTM):
    print("  NOT overwriting existing " + OUTM)
else:
    with open(OUTM, "w") as fh:
        fh.write("genome_id\tset\tgenus\tarm\tclade\ttotal\t" + "\t".join(fams) + "\n")
        for g in sorted(meta):
            r = meta[g]
            fh.write("%s\t%s\t%s\t%s\t%s\t%d\t" % (g, r["set"], r["genus"], r["arm"],
                                                   r["clade"], counts.get(g, 0)))
            fh.write("\t".join(str(gm[g].get(f, 0)) for f in fams) + "\n")
    print("  wrote %s (%d genomes x %d families)" % (OUTM, len(meta), len(fams)))

print("")
print("  dbCAN v13.0, HMMER 3.4 hmmsearch, domain i-Evalue < 1e-15 and HMM")
print("  coverage > 0.35. Search ran at E < 1e-10 so the removed tier is")
print("  visible above. Report both numbers in methods.")
print("")
print("PARSE_CAZY_FOCAL_V1_20260806 COMPLETE")
# PARSE_CAZY_FOCAL_V1_20260806
