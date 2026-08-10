# Bracken profiles from the wild samples are summarised at phylum and genus level.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/bracken_wild_bacillota.py
# Output: results/bracken_wild_phylum.tsv, results/bracken_wild_bacillota_genus.tsv
import os, csv
from collections import defaultdict

# BRACKEN_WILD_BACILLOTA V1. Replaces bracken_genus_catalog_wide.py V1, which
# had three errors: no wild-sample filter (171 samples pooled instead of 44),
# no phylum restriction on the Bracken side (so 4,317 "reads_only" genera were
# mostly not Bacillota at all), and a print that labelled the genus count as
# an SGB count. The catalog side has NO family filter, which was the original
# bug in bracken_recovery.py line 49.

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
BRACK = "/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results_bracken"
SGB = B + "/data/sgb_manifest.tsv"
MAGS = B + "/data/herptile_bacillota_A_HQ_manifest_with_source.tsv"
DB = "pluspf_20251015"
OUT_G = B + "/results/bracken_wild_bacillota_genus.tsv"
OUT_P = B + "/results/bracken_wild_phylum.tsv"

for p in (OUT_G, OUT_P):
    if os.path.exists(p):
        raise SystemExit("REFUSING TO OVERWRITE: " + p)

def read_tsv(path):
    with open(path, newline="") as fh:
        r = [x for x in csv.reader(fh, delimiter="\t") if x]
    return r[0], r[1:]

# ------------------------------------------------- WILD SAMPLES
hdr, rows = read_tsv(MAGS)
I = {c: i for i, c in enumerate(hdr)}
wild_samples = set()
all_samples = set()
for p in rows:
    s = p[I["sample_id_full"]].strip()
    all_samples.add(s)
    if p[I["source"]].strip().upper() == "WILD":
        wild_samples.add(s)

print("SAMPLE SIDE")
print("  sample_id_full in MAG manifest: %d" % len(all_samples))
print("  WILD sample_id_full:            %d" % len(wild_samples))

def gpath(s, level):
    return os.path.join(BRACK, s, s + "." + DB + ".bracken." + level + ".tsv")

wild_with_bracken = sorted(s for s in wild_samples if os.path.exists(gpath(s, "G")))
print("  WILD samples with a %s genus file: %d" % (DB, len(wild_with_bracken)))
print("  (CHATINDEX records 44 wild sample-runs with Bracken)")
missing = sorted(wild_samples - set(wild_with_bracken))
if missing:
    print("  WILD SAMPLES WITH NO BRACKEN FILE: %d" % len(missing))
    for s in missing[:10]:
        print("    " + s)
print()

def read_bracken(s, level):
    p = gpath(s, level)
    if not os.path.exists(p):
        return None
    with open(p, newline="") as fh:
        r = [x for x in csv.reader(fh, delimiter="\t") if x]
    if not r:
        return None
    h = r[0]
    try:
        ni = h.index("name")
        fi = h.index("fraction_total_reads")
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

# ------------------------------------------------- COMMUNITY CONTEXT
NS = len(wild_with_bracken)
if NS == 0:
    raise SystemExit("no wild samples with Bracken; check sample_id_full matching")

psum = defaultdict(float)
pn = defaultdict(int)
for s in wild_with_bracken:
    d = read_bracken(s, "P")
    if not d:
        continue
    for k, v in d.items():
        psum[k] += v
        pn[k] += 1

print("=" * 78)
print("COMMUNITY CONTEXT, %d WILD SAMPLES. Read every gap against this." % NS)
print("=" * 78)
print("  %-28s %10s %10s" % ("phylum", "mean_frac", "n_samples"))
for k in sorted(psum, key=lambda x: -psum[x])[:12]:
    print("  %-28s %10.4f %10d" % (k, psum[k] / NS, pn[k]))
bacillota = psum.get("Bacillota", 0.0) / NS
print()
print("  Bacillota mean read fraction: %.4f" % bacillota)
print("  Bracken reports NCBI Bacillota, which is broader than GTDB")
print("  Bacillota_A. This is an UPPER BOUND on the clade we assembled.")
print()

pf = open(OUT_P, "w")
pf.write("phylum\tmean_fraction_wild\tn_wild_samples_detected\n")
for k in sorted(psum, key=lambda x: -psum[x]):
    pf.write("%s\t%.8f\t%d\n" % (k, psum[k] / NS, pn[k]))
pf.close()
print("  wrote %s" % OUT_P)
print()

# ------------------------------------------------- CATALOG, NO FAMILY FILTER
hdr, rows = read_tsv(SGB)
J = {c: i for i, c in enumerate(hdr)}
wild_genus_family = defaultdict(set)
all_genus = set()
for p in rows:
    g = p[J["genus"]].strip()
    if not g or g.upper() == "UNASSIGNED":
        continue
    all_genus.add(g)
    if p[J["has_wild"]].strip() == "yes":
        wild_genus_family[g].add(p[J["family"]].strip())

def base(g):
    parts = g.rsplit("_", 1)
    if len(parts) == 2 and parts[1].isalpha() and parts[1].isupper() \
       and len(parts[1]) <= 2:
        return parts[0]
    return g

def is_placeholder(g):
    b = base(g)
    return not (b and b[0].isupper() and b[1:].islower() and b.isalpha())

base_to_gtdb = defaultdict(set)
for g in all_genus:
    base_to_gtdb[base(g)].add(g)

placeholder_wild = sorted(g for g in wild_genus_family if is_placeholder(g))
named_wild = sorted(g for g in wild_genus_family if not is_placeholder(g))

print("CATALOG SIDE, ALL %d SGBs, NO FAMILY FILTER" % len(rows))
print("  named genera anywhere in catalog:  %d" % len(all_genus))
print("  genera with a wild SGB:            %d" % len(wild_genus_family))
print("    testable by NCBI name matching:  %d" % len(named_wild))
print("    GTDB placeholders, UNTESTABLE:   %d" % len(placeholder_wild))
print()

# ------------------------------------------------- BACILLOTA-ONLY GENUS SET
bac_genera = set()
for s in wild_with_bracken:
    d = read_bracken(s, "F")
    if d:
        pass
gsum = defaultdict(float)
gmax = defaultdict(float)
gn = defaultdict(int)
for s in wild_with_bracken:
    d = read_bracken(s, "G")
    if not d:
        continue
    for k, v in d.items():
        gsum[k] += v
        gn[k] += 1
        if v > gmax[k]:
            gmax[k] = v

# Restrict the reads side to genera that map into the catalog's phylum by
# name. A Bracken genus is in scope only if its base name matches a GTDB
# genus present ANYWHERE in this Bacillota_A catalog. That is the only
# phylum test available without an NCBI taxonomy dump.
in_scope = [g for g in gsum if g in base_to_gtdb]
print("READS SIDE, RESTRICTED TO THE CATALOG'S PHYLUM BY NAME")
print("  genera in the wild Bracken output:        %d" % len(gsum))
print("  of those, matching a GTDB genus in this")
print("  Bacillota_A catalog:                      %d" % len(in_scope))
print("  NOTE: genera named only in NCBI with no counterpart in this catalog")
print("  cannot be separated into 'absent from the phylum' vs 'not recovered'")
print("  by name matching alone. They are reported as OUT OF SCOPE, not as")
print("  recovery failures.")
print()

B_R = defaultdict(list)
out = []
for g in sorted(in_scope):
    hits = base_to_gtdb[g]
    wild_hits = sorted(x for x in hits if x in wild_genus_family)
    fams = sorted({f for x in wild_hits for f in wild_genus_family[x]})
    if len(hits) > 1:
        b = "ambiguous_gtdb_split"
    elif wild_hits:
        b = "reads_and_wild_genome"
    else:
        b = "reads_and_captive_genome_only"
    B_R[b].append(g)
    out.append((g, gn[g], gsum[g] / NS, gmax[g], b,
                ";".join(sorted(hits)), ";".join(fams)))

print("=" * 78)
print("BUCKETS, IN-SCOPE GENERA ONLY")
print("=" * 78)
for b in ("reads_and_wild_genome", "reads_and_captive_genome_only",
          "ambiguous_gtdb_split"):
    print("  %-32s %4d" % (b, len(B_R[b])))
print("  %-32s %4d  (no NCBI name, never assignable)"
      % ("wild_genome_no_reads_possible", len(placeholder_wild)))
print()

nr = sorted([x for x in out if x[4] == "reads_and_captive_genome_only"],
            key=lambda r: -r[2])
if nr:
    print("  IN CATALOG BUT ONLY FROM CAPTIVE ANIMALS:")
    print("  %-26s %8s %12s %12s" % ("genus", "samples", "mean_frac", "max"))
    for r in nr[:15]:
        print("  %-26s %8d %12.6f %12.6f" % (r[0], r[1], r[2], r[3]))
    print()

rw = sorted([x for x in out if x[4] == "reads_and_wild_genome"],
            key=lambda r: -r[2])
print("  RECOVERED FROM WILD ANIMALS, BY READ ABUNDANCE:")
print("  %-26s %8s %12s %12s  %s" % ("genus", "samples", "mean_frac", "max",
                                     "GTDB family"))
for r in rw[:20]:
    print("  %-26s %8d %12.6f %12.6f  %s" % (r[0], r[1], r[2], r[3], r[6]))
print()

if rw:
    tot = sum(r[2] for r in rw)
    print("  summed mean read fraction of recovered genera: %.4f" % tot)
    if bacillota:
        print("  as a share of Bacillota reads: %.1f%%" % (100.0 * tot / bacillota))
print()

f = open(OUT_G, "w")
f.write("genus_ncbi\tn_wild_samples_detected\tmean_fraction_wild\t"
        "max_fraction_wild\tbucket\tgtdb_genera_matched\t"
        "gtdb_families_of_wild_sgbs\n")
for r in sorted(out, key=lambda x: -x[2]):
    f.write("%s\t%d\t%.8f\t%.8f\t%s\t%s\t%s\n" % r)
f.close()
print("  wrote %s" % OUT_G)

print()
print("BRACKEN_WILD_BACILLOTA_V1_20260806_COMPLETE")
# BRACKEN_WILD_BACILLOTA_V1_20260806_COMPLETE
