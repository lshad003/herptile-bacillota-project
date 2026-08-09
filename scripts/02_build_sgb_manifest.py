#!/usr/bin/env python3
# Builds the SGB manifest from the dRep cluster tables.
# Source: ruminococcaceae-agent/scripts/build_sgb_manifest.py
# Output: data/sgb_manifest.tsv

import os, sys
from collections import defaultdict, Counter

B = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
DREP = B + "/results/drep_herptile_95ani_2229/data_tables"
CDB = DREP + "/Cdb.csv"
WDB = DREP + "/Wdb.csv"
MAN = B + "/data/herptile_bacillota_A_HQ_manifest_with_source.tsv"
OUT = B + "/data/sgb_manifest.tsv"

# Expected counts are CHATINDEX R1 values as of 2026-08-03.
EXPECT = (("total SGBs", 1171),
          ("Ruminococcaceae SGBs", 274),
          ("Ruminococcaceae wild", 220))

for p in (CDB, MAN):
    if not os.path.exists(p):
        print("MISSING:", p); sys.exit(1)

print("SCRIPT VERSION: BUILD_SGB_MANIFEST_V2_20260803")
print()
print("dRep data_tables:")
for e in sorted(os.scandir(DREP), key=lambda x: x.name):
    print("   %-24s %d bytes" % (e.name, e.stat().st_size))


def strip_ext(x):
    for e in (".fa", ".fna", ".fasta"):
        if x.endswith(e):
            return x[:-len(e)]
    return x


mag = {}
with open(MAN) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    # 'order' is the BACTERIAL order. 'clade_order' is the HOST clade.
    # Reading 'order' into a column named host_orders was the bug fixed here.
    need = ("bin_id", "taxonomy", "source", "completeness", "contamination",
            "host_taxon", "has_metadata", "diet", "animal_type", "order",
            "clade_order")
    miss = [k for k in need if k not in h]
    if miss:
        print("MANIFEST MISSING COLUMNS: %s" % ", ".join(miss)); sys.exit(1)
    I = {k: h.index(k) for k in need}
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(I.values()) or not p[I["bin_id"]]:
            continue
        try:
            c = float(p[I["completeness"]]); x = float(p[I["contamination"]])
        except ValueError:
            continue
        t = p[I["taxonomy"]]

        def f(pre):
            for z in t.split(";"):
                z = z.strip()
                if z.startswith(pre):
                    return z[len(pre):]
            return ""

        mag[p[I["bin_id"]]] = dict(
            fam=f("f__"), gen=f("g__") or "UNASSIGNED", sp=f("s__"),
            src=p[I["source"]].strip().upper(), comp=c, cont=x,
            host=p[I["host_taxon"]].strip(), diet=p[I["diet"]].strip(),
            atype=p[I["animal_type"]].strip(),
            bord=p[I["order"]].strip(),
            hclade=p[I["clade_order"]].strip(),
            meta=p[I["has_metadata"]].strip().lower() in ("true", "yes", "1"))

n_meta = sum(1 for v in mag.values() if v["meta"])
print()
print("MAGs in manifest: %d | with metadata: %d | WITHOUT: %d"
      % (len(mag), n_meta, len(mag) - n_meta))
print("  MAGs without metadata are NOT excluded. They enter SGBs and carry")
print("  blank host fields. This is the source of the 276/274 discrepancy.")

blank_clade = sum(1 for v in mag.values() if not v["hclade"])
blank_atype = sum(1 for v in mag.values() if not v["atype"])
print("  MAGs with blank clade_order: %d | blank animal_type: %d"
      % (blank_clade, blank_atype))

print()
print("SOURCE COLUMN SEPARATION CHECK:")
bo = Counter(v["bord"] for v in mag.values() if v["bord"])
hc = Counter(v["hclade"] for v in mag.values() if v["hclade"])
print("  'order' distinct values      : %d  %s"
      % (len(bo), ", ".join(sorted(bo)[:6])))
print("  'clade_order' distinct values: %d  %s"
      % (len(hc), ", ".join(sorted(hc))))
if set(bo) & set(hc):
    print("  WARNING: the two columns share values: %s" % (set(bo) & set(hc)))
else:
    print("  Disjoint value sets, so the two columns are genuinely different.")

clu = {}
with open(CDB) as fh:
    ch = fh.readline().rstrip("\n").replace('"', "").split(",")
    cg, cs = ch.index("genome"), ch.index("secondary_cluster")
    for line in fh:
        p = line.rstrip("\n").replace('"', "").split(",")
        if len(p) > max(cg, cs):
            clu[strip_ext(p[cg])] = p[cs]
print()
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
    print("NO Wdb.csv. Representatives will be chosen by completeness minus")
    print("5x contamination, which is NOT the dRep scoring function.")

members = defaultdict(list)
for g, c in clu.items():
    if g in mag:
        members[c].append(g)
orphan = [g for g in clu if g not in mag]
print("clusters with >=1 manifest genome: %d | Cdb genomes not in manifest: %d"
      % (len(members), len(orphan)))

rows = []
nofallback = 0
for c, mem in sorted(members.items()):
    rep = winner.get(c)
    if rep not in mag:
        rep = None
    if rep is None:
        nofallback += 1
        rep = max(mem, key=lambda g: mag[g]["comp"] - 5.0 * mag[g]["cont"])
    r = mag[rep]
    srcs = sorted(set(mag[g]["src"] for g in mem))
    hosts = sorted(set(mag[g]["host"] for g in mem if mag[g]["host"]))
    diets = sorted(set(mag[g]["diet"] for g in mem if mag[g]["diet"]))
    bords = sorted(set(mag[g]["bord"] for g in mem if mag[g]["bord"]))
    hclades = sorted(set(mag[g]["hclade"] for g in mem if mag[g]["hclade"]))
    atypes = sorted(set(mag[g]["atype"] for g in mem if mag[g]["atype"]))
    rows.append(dict(sgb=c, rep=rep, n=len(mem), fam=r["fam"], gen=r["gen"],
                     sp=r["sp"], comp=r["comp"], cont=r["cont"],
                     src=",".join(srcs), wild=("WILD" in srcs),
                     hosts=";".join(hosts), nhost=len(hosts),
                     diets=";".join(diets),
                     bords=";".join(bords),
                     hclades=";".join(hclades),
                     atypes=";".join(atypes),
                     meta=all(mag[g]["meta"] for g in mem),
                     anymeta=any(mag[g]["meta"] for g in mem)))
if nofallback:
    print("representatives chosen by fallback: %d" % nofallback)

with open(OUT, "w") as f:
    f.write("sgb\trepresentative\tn_mags\tfamily\tgenus\tspecies\t"
            "rep_completeness\trep_contamination\tsources\thas_wild\t"
            "host_species\tn_host_species\tdiets\tbacterial_orders\t"
            "all_have_metadata\tany_has_metadata\thost_clades\tanimal_types\n")
    for r in rows:
        f.write("%s\t%s\t%d\t%s\t%s\t%s\t%.2f\t%.2f\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n"
                % (r["sgb"], r["rep"], r["n"], r["fam"], r["gen"], r["sp"],
                   r["comp"], r["cont"], r["src"], "yes" if r["wild"] else "no",
                   r["hosts"], r["nhost"], r["diets"], r["bords"],
                   "yes" if r["meta"] else "no",
                   "yes" if r["anymeta"] else "no",
                   r["hclades"], r["atypes"]))

print()
print("=== SGB MANIFEST ===")
print("  SGBs written              : %d" % len(rows))
print("  with >=1 wild MAG         : %d" % sum(1 for r in rows if r["wild"]))
print("  all members have metadata : %d" % sum(1 for r in rows if r["meta"]))
print("  singleton SGBs            : %d" % sum(1 for r in rows if r["n"] == 1))
print("  largest SGB               : %d MAGs" % max(r["n"] for r in rows))
print("  SGBs with blank host_clades: %d" % sum(1 for r in rows if not r["hclades"]))
print("  SGBs spanning >1 host clade: %d" % sum(1 for r in rows if ";" in r["hclades"]))

print()
print("  COLUMN CHANGE, 2026-08-03:")
print("    host_orders -> bacterial_orders. Same source column ('order' in the")
print("    herptile manifest), which holds BACTERIAL orders. It was never host")
print("    data. Only this script wrote it and nothing consumed it.")
print("    NEW host_clades, from 'clade_order', which is the actual host clade.")
print("    NEW animal_types, from 'animal_type'.")
print("  Both new columns are appended at the END, so existing column positions")
print("  are unchanged.")

print()
print("  host_clades value counts (wild SGBs):")
wildrows = [r for r in rows if r["wild"]]
for v, n in Counter(r["hclades"] for r in wildrows).most_common():
    print("     %-28s %4d" % (v if v else "<blank>", n))

print()
print("  by family (top 10):")
for fam, c in Counter(r["fam"] for r in rows).most_common(10):
    w = sum(1 for r in rows if r["fam"] == fam and r["wild"])
    print("     %-26s %4d SGBs (%d with wild)" % (fam, c, w))

print()
print("  CHECK against CHATINDEX R1 (2026-08-03):")
got = {"total SGBs": len(rows),
       "Ruminococcaceae SGBs": sum(1 for r in rows if r["fam"] == "Ruminococcaceae"),
       "Ruminococcaceae wild": sum(1 for r in rows if r["fam"] == "Ruminococcaceae"
                                   and r["wild"])}
allok = True
for lab, exp in EXPECT:
    ok = got[lab] == exp
    allok = allok and ok
    print("     %-24s got %5d  expected %5d  %s"
          % (lab, got[lab], exp, "OK" if ok else "CHECK"))
if not allok:
    print("     A CHECK here means the manifest changed. Do not overwrite the")
    print("     old file's numbers in CHATINDEX without finding out why.")

print()
print("wrote", OUT)
print()
print("EVERY downstream analysis should read this file and use")
print("`representative` as the genome, not the raw manifest.")
print("host_species is SEMICOLON delimited. Splitting on comma silently drops")
print("every multi-host SGB.")
print("DONE_SGB_MANIFEST_V2")
# SENTINEL_END
