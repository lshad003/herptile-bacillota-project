# Read classifications are joined to genome recovery across the reciprocal genera, intersected with the catalog.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/bracken_reciprocal_join_v2.py
# Output: results/bracken_reciprocal_recovery_v2.tsv
import os, csv
from collections import defaultdict, Counter

# BRACKEN_RECIPROCAL_JOIN V2. V1 applied the reciprocal filter to all 2,773
# GTDB genera in the mapping file instead of intersecting with the catalog,
# so 1,354 genera not in the catalog at all (Bacteroides, Streptomyces,
# Bradyrhizobium) landed in the "captive only" bucket with 0 SGBs, and 1,405
# unrelated genera landed in "no reads". V1's recovered table and reverse
# block were correct; the two other buckets were not.
# V2 restricts to genera with at least one SGB in the catalog, then splits on
# whether any of those SGBs come from a wild animal.

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
BRACK = "/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results_bracken"
RECIP = B + "/results/gtdb_ncbi_reciprocal.tsv"
SGB = B + "/data/sgb_manifest.tsv"
MAGS = B + "/data/herptile_bacillota_A_HQ_manifest_with_source.tsv"
DB = "pluspf_20251015"
OUT = B + "/results/bracken_reciprocal_recovery_v2.tsv"
THR_COL = "verdict_t70"

if os.path.exists(OUT):
    raise SystemExit("REFUSING TO OVERWRITE: " + OUT)

def read_tsv(path):
    with open(path, newline="") as fh:
        r = [x for x in csv.reader(fh, delimiter="\t") if x]
    return r[0], r[1:]

hdr, rows = read_tsv(RECIP)
I = {c: i for i, c in enumerate(hdr)}
to_ncbi, verdict = {}, {}
for p in rows:
    g = p[I["gtdb_genus"]].strip()
    verdict[g] = p[I[THR_COL]].strip()
    if verdict[g] == "reciprocal":
        to_ncbi[g] = p[I["best_ncbi_genus"]].strip()

hdr, rows = read_tsv(SGB)
J = {c: i for i, c in enumerate(hdr)}
wild_n, cap_n = Counter(), Counter()
fam = defaultdict(set)
for p in rows:
    g = p[J["genus"]].strip()
    if not g or g.upper() == "UNASSIGNED":
        continue
    fam[g].add(p[J["family"]].strip())
    if p[J["has_wild"]].strip() == "yes":
        wild_n[g] += 1
    else:
        cap_n[g] += 1

in_catalog = set(wild_n) | set(cap_n)
wild = set(wild_n)
captive_only = set(cap_n) - wild

print("CATALOG SIDE, %d SGBs, NO FAMILY FILTER" % len(rows))
print("  named genera in the catalog:  %d" % len(in_catalog))
print("    with a wild SGB:            %d" % len(wild))
print("    captive-derived SGBs only:  %d" % len(captive_only))
print()
print("  THE FILTER V1 LACKED: only these %d genera enter the join."
      % len(in_catalog))
print()

hdr, rows = read_tsv(MAGS)
K = {c: i for i, c in enumerate(hdr)}
wild_s = set(p[K["sample_id_full"]].strip() for p in rows
             if p[K["source"]].strip().upper() == "WILD")

def gp(s, lv):
    return os.path.join(BRACK, s, s + "." + DB + ".bracken." + lv + ".tsv")

samples = sorted(s for s in wild_s if os.path.exists(gp(s, "G")))
NS = len(samples)
print("SAMPLE SIDE: %d WILD sample_id_full, %d with a %s genus file"
      % (len(wild_s), NS, DB))
if NS == 0:
    raise SystemExit("no wild samples")

def read_b(s, lv):
    p = gp(s, lv)
    if not os.path.exists(p):
        return None
    with open(p, newline="") as fh:
        r = [x for x in csv.reader(fh, delimiter="\t") if x]
    if not r:
        return None
    h = r[0]
    try:
        ni, fi = h.index("name"), h.index("fraction_total_reads")
    except ValueError:
        return None
    d = {}
    for row in r[1:]:
        if ni < len(row) and fi < len(row):
            try:
                d[row[ni].strip()] = float(row[fi])
            except ValueError:
                pass
    return d

psum = Counter()
for s in samples:
    d = read_b(s, "P")
    if d:
        for k, v in d.items():
            psum[k] += v
bac = psum.get("Bacillota", 0.0) / NS

gsum, gn, gmax = Counter(), Counter(), defaultdict(float)
for s in samples:
    d = read_b(s, "G")
    if not d:
        continue
    for k, v in d.items():
        gsum[k] += v
        gn[k] += 1
        if v > gmax[k]:
            gmax[k] = v

print()
print("=" * 78)
print("COMMUNITY CONTEXT, %d WILD SAMPLES" % NS)
print("=" * 78)
for k, v in psum.most_common(6):
    print("  %-26s %.4f" % (k, v / NS))
print()
print("  NCBI Bacillota (%.4f) is broader than GTDB Bacillota_A and is an"
      % bac)
print("  UPPER BOUND on the clade this catalog covers.")
print()

testable = sorted(g for g in in_catalog if g in to_ncbi)
rec_w, rec_c, no_reads_w, no_reads_c = [], [], [], []
for g in testable:
    nc = to_ncbi[g]
    r = (g, nc, gn.get(nc, 0), gsum.get(nc, 0.0) / NS, gmax.get(nc, 0.0),
         ";".join(sorted(fam[g])), wild_n.get(g, 0), cap_n.get(g, 0))
    if gn.get(nc, 0) > 0:
        (rec_w if g in wild else rec_c).append(r)
    else:
        (no_reads_w if g in wild else no_reads_c).append(r)

print("=" * 78)
print("BUCKETS. %d of %d catalog genera are reciprocally testable."
      % (len(testable), len(in_catalog)))
print("=" * 78)
print("  reads detected, wild genome recovered:        %d" % len(rec_w))
print("  reads detected, genome from CAPTIVE hosts only: %d" % len(rec_c))
print("  no reads assigned, wild genome recovered:     %d" % len(no_reads_w))
print("  no reads assigned, captive genome only:       %d" % len(no_reads_c))
print()

def show(t, recs, n=40):
    if not recs:
        return
    print("  " + t)
    print("  %-24s %-22s %6s %11s %11s %6s %6s"
          % ("GTDB genus", "NCBI name", "samp", "mean_frac", "max_frac",
             "wildS", "capS"))
    for r in sorted(recs, key=lambda x: -x[3])[:n]:
        print("  %-24s %-22s %6d %11.6f %11.6f %6d %6d"
              % (r[0], r[1], r[2], r[3], r[4], r[6], r[7]))
    if len(recs) > n:
        print("    ... %d more" % (len(recs) - n))
    print()

show("READS AND A WILD GENOME:", rec_w)
show("READS, BUT EVERY GENOME CAME FROM A CAPTIVE ANIMAL:", rec_c)
show("A WILD GENOME BUT NO READS ASSIGNED:", no_reads_w)
show("CAPTIVE GENOME ONLY, NO READS ASSIGNED:", no_reads_c)

tw = sum(r[3] for r in rec_w)
tc = sum(r[3] for r in rec_c)
print("  summed mean read fraction:")
print("    genera with a wild genome:      %.4f" % tw)
print("    genera with captive genomes only: %.4f" % tc)
if bac:
    print("    wild-recovered as a share of Bacillota reads: %.1f%%"
          % (100.0 * tw / bac))
print()

print("=" * 78)
print("WHY MOST WILD CATALOG GENERA CANNOT BE TESTED AT ALL")
print("=" * 78)
vc = Counter()
for g in wild:
    v = verdict.get(g, "no_ncbi_counterpart")
    if v in ("not_in_gtdb", "not_in_reciprocal_file"):
        v = "no_ncbi_counterpart"
    vc[v] += 1
print("  wild catalog genera: %d" % len(wild))
for k in ("reciprocal", "ncbi_name_is_shared", "both_directions_fail",
          "forward_ambiguous", "no_ncbi_counterpart"):
    if vc[k]:
        print("    %-24s %4d" % (k, vc[k]))
print()
print("  NOTE: no_ncbi_counterpart means no GTDB reference genome of that")
print("  genus carries an NCBI genus name. It does NOT mean the genus is")
print("  absent from GTDB. V1 mislabelled this as not_in_gtdb.")
print()
print("  DO NOT collapse the untestable categories into 'not detected'.")
print()

f = open(OUT, "w")
f.write("gtdb_genus\tncbi_genus\tbucket\tn_wild_samples_detected\t"
        "mean_fraction_wild\tmax_fraction_wild\tgtdb_families\t"
        "wild_sgbs\tcaptive_sgbs\n")
for tag, recs in (("reads_and_wild_genome", rec_w),
                  ("reads_and_captive_genome_only", rec_c),
                  ("wild_genome_no_reads", no_reads_w),
                  ("captive_genome_no_reads", no_reads_c)):
    for r in sorted(recs, key=lambda x: -x[3]):
        f.write("%s\t%s\t%s\t%d\t%.8f\t%.8f\t%s\t%d\t%d\n"
                % (r[0], r[1], tag, r[2], r[3], r[4], r[5], r[6], r[7]))
f.close()
print("  wrote %s" % OUT)
print()
print("BRACKEN_RECIPROCAL_JOIN_V2_20260806_COMPLETE")
# BRACKEN_RECIPROCAL_JOIN_V2_20260806_COMPLETE
