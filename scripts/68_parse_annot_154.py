#!/usr/bin/env python3
# Annotation output is parsed into presence matrices and a per-genome threshold summary.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/parse_annot_154.py
# Output: results/annot_154_threshold_summary.tsv, annot_154_kofam_presence.tsv, annot_154_pfam_presence.tsv
import os, sys

BASE    = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
TARGETS = os.path.join(BASE, "work/cazy_focal/cazy_targets.tsv")
PFDIR   = os.path.join(BASE, "work/annot_154/pfam")
KODIR   = os.path.join(BASE, "work/annot_154/kofam")
OUT_SUM = os.path.join(BASE, "results/annot_154_threshold_summary.tsv")
OUT_KO  = os.path.join(BASE, "results/annot_154_kofam_presence.tsv")
OUT_PF  = os.path.join(BASE, "results/annot_154_pfam_presence.tsv")

CALIB_KO = {"K00940": "ndk", "K01937": "pyrG"}
CALIB_PF = {"PF00334": "NDK", "PF06418": "CTP_synth_N"}

def die(msg):
    sys.stderr.write("FATAL: " + msg + "\n")
    sys.exit(1)

for p in (OUT_SUM, OUT_KO, OUT_PF):
    if os.path.exists(p):
        die("output exists, refusing to overwrite: " + p)

if not os.path.isfile(TARGETS):
    die("missing " + TARGETS)

rows = []
with open(TARGETS) as fh:
    thead = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        line = line.rstrip("\n")
        if line.strip() == "":
            continue
        rows.append(line.split("\t"))

if len(rows) != 154:
    die("expected 154 target rows, found %d" % len(rows))
if len(thead) < 4:
    die("cazy_targets.tsv has %d header columns, expected at least 4" % len(thead))
if sorted(int(r[0]) for r in rows) != list(range(154)):
    die("target index column is not exactly 0-153")

def count_faa(path):
    n = 0
    with open(path) as fh:
        for line in fh:
            if line[0] == ">":
                n += 1
    return n

def parse_kofam(path):
    n_rows = n_above = n_nothr = n_bad = 0
    genes_hit = set(); genes_above = set()
    ko_above = set(); ko_seen = set()
    with open(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                if line.strip() != "":
                    n_bad += 1
                continue
            n_rows += 1
            genes_hit.add(f[1])
            ko_seen.add(f[2])
            if f[3].strip() in ("", "-"):
                n_nothr += 1
            if f[0].strip() == "*":
                n_above += 1
                genes_above.add(f[1])
                ko_above.add(f[2])
    if n_bad:
        die("%d malformed rows in %s" % (n_bad, path))
    return dict(rows=n_rows, above=n_above, nothr=n_nothr,
                genes_hit=len(genes_hit), genes_above=len(genes_above),
                ko_above=ko_above, ko_below_only=(ko_seen - ko_above))

def parse_pfam(path):
    n_rows = n_bad = 0
    genes = set(); fams = set()
    with open(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split(None, 22)
            if len(f) < 22:
                if line.strip() != "":
                    n_bad += 1
                continue
            n_rows += 1
            genes.add(f[0])
            acc = f[4]
            fams.add(acc.split(".")[0] if acc.startswith("PF") else f[3])
    if n_bad:
        die("%d malformed rows in %s" % (n_bad, path))
    return dict(rows=n_rows, genes=len(genes), fams=fams)

ko_sets = {}; pf_sets = {}; summary = []
for i, r in enumerate(rows):
    gid = r[1]; faa = r[2]
    kof = os.path.join(KODIR, gid + ".kofam.tsv")
    pff = os.path.join(PFDIR, gid + ".domtbl")
    for p in (faa, kof, pff):
        if not os.path.isfile(p):
            die("missing input for %s: %s" % (gid, p))
        if os.path.getsize(p) == 0:
            die("zero-byte input for %s: %s" % (gid, p))
    npro = count_faa(faa)
    if npro < 300:
        die("only %d proteins for %s" % (npro, gid))
    K = parse_kofam(kof)
    P = parse_pfam(pff)
    ko_sets[gid] = K["ko_above"]
    pf_sets[gid] = P["fams"]
    rec = dict(zip(thead, r))
    rec["n_proteins"] = npro
    rec["kofam_rows_total"] = K["rows"]
    rec["kofam_rows_above"] = K["above"]
    rec["kofam_rows_no_threshold"] = K["nothr"]
    rec["kofam_row_frac_above"] = round(K["above"] / K["rows"], 5) if K["rows"] else 0.0
    rec["kofam_genes_any_hit"] = K["genes_hit"]
    rec["kofam_genes_above"] = K["genes_above"]
    rec["kofam_annot_rate"] = round(K["genes_above"] / npro, 5)
    rec["kofam_ko_above"] = len(K["ko_above"])
    rec["kofam_ko_below_only"] = len(K["ko_below_only"])
    rec["pfam_rows"] = P["rows"]
    rec["pfam_genes"] = P["genes"]
    rec["pfam_annot_rate"] = round(P["genes"] / npro, 5)
    rec["pfam_families"] = len(P["fams"])
    summary.append(rec)
    if (i + 1) % 25 == 0:
        sys.stderr.write("parsed %d / 154\n" % (i + 1))

scols = list(thead) + ["n_proteins", "kofam_rows_total", "kofam_rows_above",
    "kofam_rows_no_threshold", "kofam_row_frac_above", "kofam_genes_any_hit",
    "kofam_genes_above", "kofam_annot_rate", "kofam_ko_above",
    "kofam_ko_below_only", "pfam_rows", "pfam_genes", "pfam_annot_rate",
    "pfam_families"]
with open(OUT_SUM, "w") as out:
    out.write("\t".join(scols) + "\n")
    for rec in summary:
        out.write("\t".join(str(rec.get(c, "")) for c in scols) + "\n")

def write_matrix(path, sets, label):
    feats = sorted(set().union(*sets.values())) if sets else []
    with open(path, "w") as out:
        out.write("genome\t" + "\t".join(feats) + "\n")
        for rec in summary:
            g = rec[thead[1]]
            s = sets[g]
            out.write(g + "\t" + "\t".join("1" if f in s else "0" for f in feats) + "\n")
    sys.stderr.write("%s: 154 genomes x %d %s\n" % (os.path.basename(path), len(feats), label))

write_matrix(OUT_KO, ko_sets, "KOs")
write_matrix(OUT_PF, pf_sets, "Pfam families")

setcol = thead[3]
groups = {}
for rec in summary:
    groups.setdefault(rec[setcol], []).append(rec)

def med(v):
    v = sorted(v); n = len(v)
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2.0

print("")
print("=== SET COLUMN: %s ===" % setcol)
for k in sorted(groups):
    print("  %-14s n=%d" % (k, len(groups[k])))
print("")
print("=== KOFAM THRESHOLD BEHAVIOUR, all 154 ===")
tot = sum(r["kofam_rows_total"] for r in summary)
abv = sum(r["kofam_rows_above"] for r in summary)
print("  rows total %d, above threshold %d (%.2f%%)" % (tot, abv, 100.0 * abv / tot))
print("  rows with no defined threshold: %d" % sum(r["kofam_rows_no_threshold"] for r in summary))
print("  median per-genome row fraction above: %.4f" % med([r["kofam_row_frac_above"] for r in summary]))
print("")
print("=== ANNOTATION RATE, fraction of proteins with a call ===")
for lab, key in (("kofam", "kofam_annot_rate"), ("pfam ", "pfam_annot_rate")):
    v = [r[key] for r in summary]
    print("  %s  median %.4f  min %.4f  max %.4f" % (lab, med(v), min(v), max(v)))
for k in sorted(groups):
    v = [r["kofam_annot_rate"] for r in groups[k]]
    w = [r["pfam_annot_rate"] for r in groups[k]]
    print("  %-14s kofam median %.4f | pfam median %.4f" % (k, med(v), med(w)))
print("")
print("=== CALIBRATION GENES, genomes carrying the call ===")
for ko, name in sorted(CALIB_KO.items()):
    print("  %s (%s), KofamScan above threshold:" % (ko, name))
    for k in sorted(groups):
        n = sum(1 for r in groups[k] if ko in ko_sets[r[thead[1]]])
        print("      %-14s %d of %d" % (k, n, len(groups[k])))
for pf, name in sorted(CALIB_PF.items()):
    print("  %s (%s), Pfam 38.2 cut_ga:" % (pf, name))
    for k in sorted(groups):
        n = sum(1 for r in groups[k] if pf in pf_sets[r[thead[1]]])
        print("      %-14s %d of %d" % (k, n, len(groups[k])))
print("")
print("ANCHORS FOR COMPARISON, from prior work on the 125 focal genomes only:")
print("  ndk   eggNOG 0.02 amphibian / 0.12 reference; PF00334 direct scan 5 of 125")
print("  pyrG  eggNOG 0/46 amph Anaerotruncus, 3/52 amph UBA866; PF06418 24 of 125;")
print("        DIAMOND 25 of 125 with positive controls 125 of 125")
print("  These anchors are 125 genomes. Do not compare them to a 154 denominator.")
print("")
print("WROTE:")
for p in (OUT_SUM, OUT_KO, OUT_PF):
    print("  " + p)
# PARSE_ANNOT_154_V1_20260808
