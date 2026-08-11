# Groups are retained only where the direction of difference agrees within both genera.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/genus_filter_merged.py
# Output: work/focal_genus_pangenome/matrices/happi_og_genus_filtered.tsv
# GENUS_FILTER_MERGED_V1_20260805
import os
from collections import Counter, defaultdict

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
WORK = os.path.join(ROOT, "work/focal_genus_pangenome")
MATD = os.path.join(WORK, "matrices")
RES = os.path.join(MATD, "happi_results_og_bacteria.tsv")
PRES = os.path.join(MATD, "presence_og_bacteria.tsv")
MET = os.path.join(MATD, "unit_metadata.tsv")
C2OG = os.path.join(MATD, "cluster_to_og.tsv")
ANN = os.path.join(WORK, "eggnog/focal.emapper.annotations")
OUT = os.path.join(MATD, "happi_og_genus_filtered.tsv")

Q = 0.05
MIN_PER_CELL = 5   # minimum genomes in a genus x group cell to judge direction


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


def tofloat(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


meta = {r["genome"]: r for r in read_tsv(MET)}
cells = defaultdict(list)
for g, r in meta.items():
    cells[(r["genus"], r["group"])].append(g)
print("=" * 74)
print("DESIGN CELLS")
print("=" * 74)
for k in sorted(cells):
    print("  %-16s %-12s n=%d" % (k[0], k[1], len(cells[k])))
GENERA = sorted({r["genus"] for r in meta.values()})

pres = {}
with open(PRES) as fh:
    gs = fh.readline().rstrip("\n").split("\t")[1:]
    for line in fh:
        f = line.rstrip("\n").split("\t")
        pres[f[0]] = frozenset(g for g, v in zip(gs, f[1:]) if v == "1")
print("  presence rows: %d" % len(pres))

res = read_tsv(RES)
sig = [r for r in res if (tofloat(r["q"]) or 1) < Q]
print("  significant at q<%.2f: %d" % (Q, len(sig)))

# ---------------------------------------------------------------- filter
print()
print("=" * 74)
print("WITHIN-GENUS CONSISTENCY")
print("=" * 74)
print("  The pooled test compares 98 amphibian against 26 reference genomes")
print("  across two genera. A group present only in the amphibian members of")
print("  ONE genus is a lineage signal, not a host signal.")
print("  Strong filter: the direction of the amphibian-minus-reference")
print("  difference must hold within EACH genus separately.")
print("  Weak filter: the group is merely present in both genera.")

rows = []
for r in sig:
    og = r["og"]
    carriers = pres.get(og, frozenset())
    diff = tofloat(r["diff"])
    per = {}
    testable = True
    for gen in GENERA:
        a = cells.get((gen, "amphibian"), [])
        b = cells.get((gen, "reference"), [])
        if len(a) < MIN_PER_CELL or len(b) < MIN_PER_CELL:
            testable = False
            per[gen] = (None, None, None)
            continue
        pa = sum(1 for x in a if x in carriers) / float(len(a))
        pr = sum(1 for x in b if x in carriers) / float(len(b))
        per[gen] = (pa, pr, pa - pr)
    signs = [per[g][2] for g in GENERA if per[g][2] is not None]
    consistent = (testable and len(signs) == len(GENERA)
                  and all((s > 0) == (diff > 0) for s in signs)
                  and all(abs(s) > 0 for s in signs))
    amph_gen = {meta[x]["genus"] for x in carriers if meta[x]["group"] == "amphibian"}
    rows.append(dict(r=r, og=og, diff=diff, per=per,
                     consistent=consistent, both_present=len(amph_gen) >= 2,
                     carriers=carriers))

amph_sig = [x for x in rows if x["diff"] > 0]
ref_sig = [x for x in rows if x["diff"] < 0]
print()
print("  %-22s %8s %14s %16s" % ("direction", "n", "both genera", "direction holds"))
for lab, sub in (("amphibian-higher", amph_sig), ("reference-higher", ref_sig)):
    print("  %-22s %8d %14d %16d"
          % (lab, len(sub), sum(1 for x in sub if x["both_present"]),
             sum(1 for x in sub if x["consistent"])))

# ---------------------------------------------------------------- annotate
og_clusters = defaultdict(list)
for r in read_tsv(C2OG):
    og_clusters[r["og_bacteria"]].append(r["cluster"])

ann, hdr = {}, None
with open(ANN) as fh:
    for line in fh:
        if line.startswith("#query"):
            hdr = line.lstrip("#").rstrip("\n").split("\t")
            continue
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if hdr and len(f) >= len(hdr):
            ann[f[0]] = dict(zip(hdr, f))


def og_ann(og):
    names, descs, cogs, kos, cazy = Counter(), Counter(), Counter(), Counter(), Counter()
    for cid in og_clusters.get(og, []):
        a = ann.get(cid)
        if not a:
            continue
        for src, tgt in ((a.get("Preferred_name", ""), names),
                         (a.get("Description", ""), descs),
                         (a.get("COG_category", ""), cogs),
                         (a.get("KEGG_ko", ""), kos),
                         (a.get("CAZy", ""), cazy)):
            v = src.strip()
            if v and v != "-":
                tgt[v] += 1
    pick = lambda c: c.most_common(1)[0][0] if c else "-"
    return pick(names), pick(descs), pick(cogs), pick(kos), pick(cazy)


keep = [x for x in rows if x["consistent"]]
keep.sort(key=lambda x: tofloat(x["r"]["p"]))
print()
print("=" * 74)
print("GROUPS PASSING THE WITHIN-GENUS DIRECTION FILTER: %d" % len(keep))
print("=" * 74)

with open(OUT, "w") as f:
    f.write("og\tdirection\tprev_amphibian\tprev_reference\tdiff\tp\tq\t"
            + "\t".join("%s_amph\t%s_ref\t%s_diff" % (g, g, g) for g in GENERA)
            + "\tboth_genera\tdirection_consistent\tgene\tdescription\tCOG\tKEGG_ko\tCAZy\n")
    for x in rows:
        r = x["r"]
        nm, ds, cg, ko, cz = og_ann(x["og"])
        cellstr = "\t".join(
            ("%.3f\t%.3f\t%.3f" % x["per"][g]) if x["per"][g][0] is not None
            else "NA\tNA\tNA" for g in GENERA)
        f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                % (x["og"], "amphibian" if x["diff"] > 0 else "reference",
                   r["prev_amphibian"], r["prev_reference"], r["diff"],
                   r["p"], r["q"], cellstr,
                   "yes" if x["both_present"] else "no",
                   "yes" if x["consistent"] else "no",
                   nm, ds[:80], cg, ko[:40], cz))
print("  wrote %s" % OUT)

for direction in ("amphibian", "reference"):
    sub = [x for x in keep if (x["diff"] > 0) == (direction == "amphibian")]
    print()
    print("  %s-HIGHER, top 12 by p (%d total)" % (direction.upper(), len(sub)))
    print("    %-12s %6s %6s %8s  %-28s %s"
          % ("og", "amph", "ref", "q", "gene / description", "per-genus diff"))
    for x in sub[:12]:
        r = x["r"]
        nm, ds, cg, ko, cz = og_ann(x["og"])
        lab = nm if nm != "-" else (ds[:26] if ds != "-" else "unannotated")
        pg = " ".join("%s%+.2f" % (g[:4], x["per"][g][2]) for g in GENERA)
        print("    %-12s %6.2f %6.2f %8.2g  %-28s %s"
              % (x["og"][:12], tofloat(r["prev_amphibian"]),
                 tofloat(r["prev_reference"]), tofloat(r["q"]), lab[:28], pg))

print()
print("=" * 74)
print("SUMMARY BY COG CATEGORY, filtered set")
print("=" * 74)
for direction in ("amphibian", "reference"):
    sub = [x for x in keep if (x["diff"] > 0) == (direction == "amphibian")]
    cc = Counter()
    for x in sub:
        _, _, cg, _, _ = og_ann(x["og"])
        for ch in cg:
            if ch.isalpha():
                cc[ch] += 1
    print("  %-10s %s" % (direction, dict(cc.most_common(10)) if cc else "none"))
print()
print("  COG letters: C energy, E amino acid, F nucleotide, G carbohydrate,")
print("  H coenzyme, J translation, K transcription, L replication,")
print("  M cell wall, P inorganic ion, S unknown, T signalling, V defence.")

# SENTINEL_END
