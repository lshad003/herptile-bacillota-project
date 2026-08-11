# GTDB r220 reference alignments are extracted and masked for the joint tree.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/extract_ref_msa_rum.py
# Output: work/rep_tree/gtdb_ref_rum_masked.faa, gtdb_ref_rum_metadata.tsv
# EXTRACT_REF_MSA_RUM_V1_20260804
import os, sys
from collections import Counter

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
REF_MSA = "/srv/projects/db/gtdbtk/220/msa/gtdb_r220_bac120.faa"
MASK = "/srv/projects/db/gtdbtk/220/masks/gtdb_r220_bac120.mask"
REFS = os.path.join(ROOT, "data/tasks/gtdb_ruminococcaceae_paths.tsv")
GTDB_TAX = "/srv/projects/db/gtdbtk/220/taxonomy/gtdb_taxonomy.tsv"

OUTD = os.path.join(ROOT, "work/rep_tree")
OUT_FA = os.path.join(OUTD, "gtdb_ref_rum_masked.faa")
OUT_META = os.path.join(OUTD, "gtdb_ref_rum_metadata.tsv")

EXPECT_UNMASKED = 41084
EXPECT_MASKED = 5035

if not os.path.isdir(OUTD):
    os.makedirs(OUTD)

# ------------------------------------------------------------------ mask
with open(MASK) as fh:
    mask = fh.readline().strip()
keep = [i for i, c in enumerate(mask) if c == "1"]
print("=" * 74)
print("MASK")
print("=" * 74)
print("  length: %d (expect %d)" % (len(mask), EXPECT_UNMASKED))
print("  ones  : %d (expect %d)" % (len(keep), EXPECT_MASKED))
print("  characters present: %s" % dict(Counter(mask)))
if len(mask) != EXPECT_UNMASKED or len(keep) != EXPECT_MASKED:
    sys.exit("MASK DOES NOT MATCH EXPECTED DIMENSIONS. Stopping.")
print("  Applying this mask reproduces exactly what GTDB-Tk's align step did")
print("  for the user MSAs, so masked reference sequences land in the same")
print("  5,035-column space and are concatenable with them.")

# ------------------------------------------------------------ targets
want = set()
with open(REFS) as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f and f[0].strip():
            want.add(f[0].strip())
print()
print("  Ruminococcaceae accessions wanted: %d" % len(want))

genus = {}
with open(GTDB_TAX) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        acc = f[0].strip().replace("GB_", "").replace("RS_", "")
        if acc in want:
            g = ""
            for z in f[1].split(";"):
                z = z.strip()
                if z.startswith("g__"):
                    g = z[3:]
            genus[acc] = g

# --------------------------------------------------- stream the 4.4 GB msa
print()
print("=" * 74)
print("STREAMING THE REFERENCE MSA (4.4 GB, never loaded into memory)")
print("=" * 74)
found = {}
cur_acc, cur_buf, keeping = None, [], False
n_seen = 0
bad_len = []


def flush():
    if cur_acc is None or not keeping:
        return
    s = "".join(cur_buf)
    if len(s) != EXPECT_UNMASKED:
        bad_len.append((cur_acc, len(s)))
        return
    found[cur_acc] = "".join(s[i] for i in keep)


with open(REF_MSA) as fh:
    for line in fh:
        if line.startswith(">"):
            flush()
            n_seen += 1
            raw = line[1:].strip().split()[0]
            cur_acc = raw.replace("GB_", "").replace("RS_", "")
            keeping = cur_acc in want
            cur_buf = []
        elif keeping:
            cur_buf.append(line.strip())
    flush()

print("  sequences scanned : %d" % n_seen)
print("  matched and masked: %d of %d wanted" % (len(found), len(want)))
missing = sorted(want - set(found))
print("  not found in MSA  : %d %s" % (len(missing), missing[:5]))
if bad_len:
    print("  WRONG UNMASKED LENGTH: %d %s" % (len(bad_len), bad_len[:3]))

if found:
    w = set(len(v) for v in found.values())
    print("  masked widths: %s" % w)
    if w != {EXPECT_MASKED}:
        sys.exit("masked width mismatch, stopping")
    gaps = sorted(v.count("-") / float(EXPECT_MASKED) for v in found.values())
    print("  gap fraction: min %.3f median %.3f max %.3f"
          % (gaps[0], gaps[len(gaps) // 2], gaps[-1]))

with open(OUT_FA, "w") as fa, open(OUT_META, "w") as mp:
    mp.write("tip\tarm\tgenome\thost_category\thost_detail\tfamily\tgenus\t"
             "n_in_cluster_this_arm\tgap_fraction\n")
    for acc in sorted(found):
        s = found[acc]
        tip = "gtdb_ref|%s" % acc
        fa.write(">%s\n%s\n" % (tip, s))
        mp.write("%s\tgtdb_ref\t%s\tnone\t\tRuminococcaceae\t%s\t1\t%.4f\n"
                 % (tip, acc, genus.get(acc, ""),
                    s.count("-") / float(EXPECT_MASKED)))

print()
print("  wrote %s (%d seqs)" % (OUT_FA, len(found)))
print("  wrote %s" % OUT_META)
print()
print("  These are NOT for the gate diagnostic: references have no host, so")
print("  they cannot break the host/clade confound. They are for a")
print("  reference-anchored figure tree showing where the herptile lineages")
print("  sit within the family.")

# SENTINEL_END
