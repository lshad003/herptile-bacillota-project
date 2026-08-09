#!/usr/bin/env python3
# Parsing of barrnap and tRNAscan-SE output for both arms.
#
# Source: ruminococcaceae-agent/scripts/parse_rrna_trna_both_arms.py
# Reads:  results/rrna_trna/ (1,171 SGB representatives)
#         results/rrna_trna_refs/ (1,246 GTDB r220 references)
# Writes: results/rrna_trna_catalog_per_genome_merged.tsv
#         results/rrna_trna_refs_per_genome.tsv
#
# Both arms are parsed by one code path so that detection criteria cannot
# differ between them.
#
# tRNA isotypes are counted by amino acid. Ile2 is counted as isoleucine and
# fMet and iMet as methionine, since these are isoacceptors of standard amino
# acids. SeC, suppressor and undetermined tRNAs are excluded from the
# twenty-amino-acid count. An earlier parser,
# ruminococcaceae-agent/scripts/parse_rrna_trna.py, matched isotype names
# against a whitelist of the twenty canonical three-letter abbreviations,
# which discards Ile2 and fMet and understates tRNA coverage; its output
# should not be used.
#
# The catalog-arm output block was added on 2026-08-08. The original script
# summarised both arms but wrote only the reference arm to disk.
import os, sys, csv
from collections import Counter, defaultdict

BASE = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
QRY  = os.path.join(BASE, "results/rrna_trna")
REF  = os.path.join(BASE, "results/rrna_trna_refs")
MANI = os.path.join(BASE, "data/sgb_manifest.tsv")
OUT  = os.path.join(BASE, "results/rrna_trna_refs_per_genome.tsv")
GAPS = os.path.join(BASE, "results/rrna_trna_refs_missing.txt")

EXPECT_REF = 1247
EXPECT_QRY = 1171
STD20 = set("Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val".split())

def die(m):
    print("")
    print("!" * 72)
    print("FAILED: " + m)
    print("!" * 72)
    sys.exit(1)

def gid(fn):
    s = fn
    for ext in (".gff", ".tsv", ".txt", ".out", ".trna", ".tRNA"):
        if s.endswith(ext):
            s = s[: -len(ext)]
    for pre in ("gtdbref__", "gtdb_ref__", "herp__", "amph__", "sgb__"):
        if s.startswith(pre):
            s = s[len(pre):]
    return s

def scan(d):
    if not os.path.isdir(d):
        return {}
    return {gid(e.name): e.path for e in os.scandir(d) if e.is_file()}

def parse_barrnap(path):
    """Return dict of gene -> 'complete' | 'partial'. Complete wins."""
    out = {}
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                attr = f[8]
                name = None
                for kv in attr.split(";"):
                    if kv.startswith("Name="):
                        name = kv[5:].strip()
                        break
                if not name:
                    continue
                g = name.replace("_rRNA", "").strip()
                st = "partial" if "partial" in attr.lower() else "complete"
                if out.get(g) != "complete":
                    out[g] = st
    except Exception as e:
        return None
    return out

def parse_trnascan(path):
    """Return (n_tRNA, set of standard AAs, Counter of raw types)."""
    n, aas, raw = 0, set(), Counter()
    try:
        with open(path) as fh:
            for line in fh:
                if not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 5:
                    f = line.split()
                if len(f) < 5:
                    continue
                if f[0].strip().lower().startswith(("sequence", "name", "----")):
                    continue
                try:
                    int(f[1].strip())
                except Exception:
                    continue
                t = f[4].strip()
                raw[t] += 1
                low = " ".join(f).lower()
                if "pseudo" in low:
                    continue
                n += 1
                if t == "Ile2":
                    t = "Ile"
                if t in ("fMet", "iMet"):
                    t = "Met"
                if t in STD20:
                    aas.add(t)
    except Exception:
        return None
    return n, aas, raw

def collect(root, label, expect):
    print("")
    print("=" * 72)
    print("ARM: %s   %s" % (label, root))
    print("=" * 72)
    bp = scan(os.path.join(root, "barrnap"))
    tp = scan(os.path.join(root, "trnascan"))
    print("  barrnap files : %d" % len(bp))
    print("  trnascan files: %d" % len(tp))
    print("  expected      : %d" % expect)
    if bp:
        print("  example id    : %s" % list(bp.keys())[0])

    empty_b = [g for g, p in bp.items() if os.path.getsize(p) == 0]
    empty_t = [g for g, p in tp.items() if os.path.getsize(p) == 0]
    print("  ZERO-BYTE barrnap : %d" % len(empty_b))
    print("  ZERO-BYTE trnascan: %d" % len(empty_t))

    ids = set(bp) | set(tp)
    miss_b = sorted((ids - set(bp)) | set(empty_b))
    miss_t = sorted((ids - set(tp)) | set(empty_t))
    print("  genomes seen (union): %d" % len(ids))
    print("  missing/empty barrnap : %d" % len(miss_b))
    print("  missing/empty trnascan: %d" % len(miss_t))
    if len(ids) != expect:
        print("  *** COUNT MISMATCH: %d seen against %d expected ***" % (len(ids), expect))

    rows, badb, badt = {}, 0, 0
    typeseen = Counter()
    for g in sorted(ids):
        rr = parse_barrnap(bp[g]) if g in bp and os.path.getsize(bp[g]) > 0 else None
        tr = parse_trnascan(tp[g]) if g in tp and os.path.getsize(tp[g]) > 0 else None
        if rr is None:
            badb += 1
        if tr is None:
            badt += 1
        else:
            typeseen.update(tr[2])
        rows[g] = (rr if rr is not None else {},
                   tr[0] if tr else -1,
                   tr[1] if tr else set(),
                   (rr is not None), (tr is not None))
    if badb or badt:
        print("  UNPARSEABLE barrnap %d, trnascan %d" % (badb, badt))
    print("  tRNA types observed (top 8): %s" % typeseen.most_common(8))
    nonstd = [t for t in typeseen if t not in STD20 and t not in ("Ile2", "fMet", "iMet")]
    if nonstd:
        print("  non-standard types excluded from the 20-AA count: %s" % sorted(nonstd)[:12])
    return rows, miss_b, miss_t

def summarise(rows, label, subset=None):
    if subset is not None:
        rows = {g: v for g, v in rows.items() if g in subset}
    usable = {g: v for g, v in rows.items() if v[3] and v[4]}
    n = len(usable)
    print("")
    print("-" * 72)
    print("%s   n = %d genomes with BOTH tools parsed" % (label, n))
    print("-" * 72)
    if n == 0:
        print("  nothing to summarise")
        return None
    def cnt(gene, state):
        return sum(1 for v in usable.values() if v[0].get(gene) == state)
    for gene in ("5S", "16S", "23S"):
        c, p = cnt(gene, "complete"), cnt(gene, "partial")
        print("  %-4s complete %5d (%5.1f%%)   partial-only %4d" % (gene, c, 100.0 * c / n, p))
    all3 = sum(1 for v in usable.values()
               if all(v[0].get(g) == "complete" for g in ("5S", "16S", "23S")))
    print("  ALL THREE complete       %5d (%5.1f%%)" % (all3, 100.0 * all3 / n))
    aa18 = sum(1 for v in usable.values() if len(v[2]) >= 18)
    print("  tRNAs for >=18 of 20 AAs %5d (%5.1f%%)" % (aa18, 100.0 * aa18 / n))
    counts = sorted(v[1] for v in usable.values() if v[1] >= 0)
    if counts:
        print("  tRNA count median %d (min %d, max %d)"
              % (counts[len(counts) // 2], counts[0], counts[-1]))
    dis = sorted(len(v[2]) for v in usable.values())
    print("  distinct AAs median %d" % dis[len(dis) // 2])
    both = sum(1 for v in usable.values()
               if all(v[0].get(g) == "complete" for g in ("5S", "16S", "23S"))
               and len(v[2]) >= 18)
    print("  MIMAG rRNA+tRNA BOTH     %5d (%5.1f%%)" % (both, 100.0 * both / n))
    return n, all3, aa18, both

refrows, rmb, rmt = collect(REF, "GTDB REFERENCE", EXPECT_REF)
qryrows, qmb, qmt = collect(QRY, "HERPTILE CATALOG", EXPECT_QRY)

rumino = None
if os.path.exists(MANI):
    with open(MANI) as fh:
        d = "," if fh.readline().count(",") > 1 else "\t"
    mrows = list(csv.DictReader(open(MANI), delimiter=d))
    hdr = list(mrows[0].keys())
    low = {h.lower(): h for h in hdr}
    kid = next((low[k] for k in ("representative", "genome", "genome_id") if k in low), None)
    kfa = next((low[k] for k in ("family", "gtdb_family") if k in low), None)
    if kid and kfa:
        rumino = set(r[kid].strip() for r in mrows
                     if "Ruminococcaceae" in (r[kfa] or ""))
        print("")
        print("  manifest Ruminococcaceae representatives: %d" % len(rumino))
        print("  of which present in the parsed catalog  : %d"
              % len(rumino & set(qryrows)))

print("")
print("#" * 72)
print("# SIDE BY SIDE, SAME PARSING CODE ON BOTH ARMS")
print("#" * 72)
R = summarise(refrows, "GTDB REFERENCE Ruminococcaceae (isolates + MAGs)")
Q = summarise(qryrows, "HERPTILE CATALOG, all SGB representatives")
QR = summarise(qryrows, "HERPTILE CATALOG, Ruminococcaceae only", rumino) if rumino else None

if R and Q:
    print("")
    print("=" * 72)
    print("THE COMPARISON")
    print("=" * 72)
    print("  all three rRNA complete:  reference %5.1f%%   catalog %5.1f%%"
          % (100.0 * R[1] / R[0], 100.0 * Q[1] / Q[0]))
    print("  >=18 AA tRNA coverage  :  reference %5.1f%%   catalog %5.1f%%"
          % (100.0 * R[2] / R[0], 100.0 * Q[2] / Q[0]))
    print("  MIMAG rRNA+tRNA both   :  reference %5.1f%%   catalog %5.1f%%"
          % (100.0 * R[3] / R[0], 100.0 * Q[3] / Q[0]))
    print("")
    print("  READ THIS BEFORE WRITING ANYTHING FROM IT:")
    print("  91.9%% of GTDB reference Ruminococcaceae are THEMSELVES MAGs, so a")
    print("  low reference rate is NOT an isolate-versus-MAG contrast. It is a")
    print("  short-read-assembly contrast. If you want the isolate comparison")
    print("  you must split the reference arm on provenance first, using")
    print("  results/gtdb_reference_provenance.tsv.")
    print("  Reference median N50 is 37,412 against 17,350 in the catalog, so")
    print("  the arms are NOT contiguity-matched. State that with the numbers.")

if os.path.exists(OUT):
    print("")
    print("  NOT overwriting existing " + OUT)
else:
    with open(OUT, "w") as fh:
        fh.write("genome\ts5\ts16\ts23\tall_three_complete\tn_trna\tn_distinct_aa\taa18\tbarrnap_ok\ttrnascan_ok\n")
        for g in sorted(refrows):
            rr, nt, aas, bok, tok = refrows[g]
            a3 = 1 if all(rr.get(x) == "complete" for x in ("5S", "16S", "23S")) else 0
            fh.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n"
                     % (g, rr.get("5S", "absent"), rr.get("16S", "absent"),
                        rr.get("23S", "absent"), a3, nt, len(aas),
                        1 if len(aas) >= 18 else 0, int(bok), int(tok)))
    print("")
    print("  wrote " + OUT)

OUT_QRY = os.path.join(BASE, "results/rrna_trna_catalog_per_genome_merged.tsv")
if os.path.exists(OUT_QRY):
    print("")
    print("  NOT overwriting existing " + OUT_QRY)
else:
    with open(OUT_QRY, "w") as fh:
        fh.write("genome\ts5\ts16\ts23\tall_three_complete\tn_trna\tn_distinct_aa\taa18\tbarrnap_ok\ttrnascan_ok\n")
        for g in sorted(qryrows):
            rr, nt, aas, bok, tok = qryrows[g]
            a3 = 1 if all(rr.get(x) == "complete" for x in ("5S", "16S", "23S")) else 0
            fh.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n"
                     % (g, rr.get("5S", "absent"), rr.get("16S", "absent"),
                        rr.get("23S", "absent"), a3, nt, len(aas),
                        1 if len(aas) >= 18 else 0, int(bok), int(tok)))
    print("")
    print("  wrote " + OUT_QRY)

gaps = sorted(set(rmb) | set(rmt))
if gaps:
    if os.path.exists(GAPS):
        print("  NOT overwriting existing " + GAPS)
    else:
        with open(GAPS, "w") as fh:
            for g in gaps:
                fh.write(g + "\n")
        print("  wrote %d gap ids to %s" % (len(gaps), GAPS))
else:
    print("  no gaps to rerun")

print("")
print("PARSE_RRNA_TRNA_BOTH_ARMS_V1_20260806 COMPLETE")
# PARSE_RRNA_TRNA_BOTH_ARMS_V1_20260806
