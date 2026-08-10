# GTDB-Tk quality flags on the wild representatives are audited against manifest contamination and GUNC calls.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/audit_wild_classify_flags.py
# Output: results/wild_classify_flagged_genomes.tsv
# AUDIT_WILD_FLAGS_V1_20260803
import os
from collections import Counter

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
SUM = os.path.join(ROOT, "results/gtdbtk_wild_sgb_r220/gtdbtk.bac120.summary.tsv")
MANIFEST = os.path.join(ROOT, "data/sgb_manifest.tsv")
GUNC = os.path.join(ROOT, "results/gunc_audit_by_sgb.tsv")
OUT = os.path.join(ROOT, "results/wild_classify_flagged_genomes.tsv")

MSA_FLOOR = 50.0
MULTIHIT = "markers with multiple hits"


def read_tsv(path):
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                f = f + [""] * (len(header) - len(f))
            rows.append(dict(zip(header, f)))
    return rows, header


def parse_tax(s):
    out = {}
    for part in s.strip().split(";"):
        part = part.strip()
        if len(part) > 3 and part[1:3] == "__":
            out[part[0]] = part[3:]
    return out


def tofloat(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def pct(a, b):
    return (100.0 * a / b) if b else 0.0


new, _ = read_tsv(SUM)
man_rows, _ = read_tsv(MANIFEST)
man = {r["representative"].strip(): r for r in man_rows}

gunc = {}
gkeys = []
if os.path.exists(GUNC):
    grows, gheader = read_tsv(GUNC)
    gkeys = gheader
    idcol = None
    for cand in ("representative", "genome", "sgb", "user_genome"):
        if cand in gheader:
            idcol = cand
            break
    print("GUNC id column: %s | columns: %s" % (idcol, ", ".join(gheader[:8])))
    if idcol:
        for r in grows:
            gunc[r[idcol].strip()] = r
else:
    print("GUNC file missing, cross-reference skipped")

# ------------------------------------------------------------ collect flags
flagged = {}
for r in new:
    g = r["user_genome"].strip()
    w = r["warnings"].strip()
    msa = tofloat(r["msa_percent"])
    reasons = []
    if w and w != "N/A":
        reasons.append("warning: " + w)
    if msa is not None and msa < MSA_FLOOR:
        reasons.append("msa_percent %.1f below %.0f" % (msa, MSA_FLOOR))
    if reasons:
        flagged[g] = dict(row=r, msa=msa, reasons=reasons,
                          multihit=(MULTIHIT in w),
                          ani_radius=("ANI radius" in w))

print()
print("=" * 74)
print("FLAGGED GENOMES")
print("=" * 74)
print("  total flagged: %d of %d (%.1f%%)" % (len(flagged), len(new), pct(len(flagged), len(new))))
print("    multiple-hit marker warning : %d" % sum(1 for v in flagged.values() if v["multihit"]))
print("    outside ANI radius warning  : %d" % sum(1 for v in flagged.values() if v["ani_radius"]))
print("    msa_percent below %.0f       : %d"
      % (MSA_FLOOR, sum(1 for v in flagged.values()
                        if v["msa"] is not None and v["msa"] < MSA_FLOOR)))
both = [g for g, v in flagged.items()
        if v["multihit"] and v["msa"] is not None and v["msa"] < MSA_FLOOR]
print("    BOTH multi-hit and low MSA  : %d %s" % (len(both), both if both else ""))

# ---------------------------------------------------- per-genome detail table
print()
print("=" * 74)
print("DETAIL, sorted by msa_percent ascending")
print("=" * 74)
print("  %-34s %6s %5s %6s %6s %5s %s"
      % ("genome", "msa%", "nmags", "comp", "cont", "GUNC", "family"))
out = open(OUT, "w")
out.write("genome\tsgb\tfamily\tgenus\tn_mags\tcompleteness\tcontamination\t"
          "msa_percent\tgunc_pass\tgunc_css\tgunc_rrs\tflags\n")

rows_sorted = sorted(flagged.items(),
                     key=lambda kv: (kv[1]["msa"] if kv[1]["msa"] is not None else 999))
gunc_fail = 0
for g, v in rows_sorted:
    r = v["row"]
    t = parse_tax(r["classification"])
    m = man.get(g, {})
    gu = gunc.get(g, {})
    gpass = gu.get("pass_gunc", "").strip() if gu else "no_gunc_row"
    gcss = gu.get("clade_separation_score", "").strip() if gu else ""
    grrs = gu.get("reference_representation_score", "").strip() if gu else ""
    if gpass and gpass not in ("pass", "no_gunc_row"):
        gunc_fail += 1
    print("  %-34s %6s %5s %6s %6s %5s %s"
          % (g[:34],
             ("%.1f" % v["msa"]) if v["msa"] is not None else "?",
             m.get("n_mags", "?"),
             m.get("rep_completeness", "?"),
             m.get("rep_contamination", "?"),
             gpass[:5],
             t.get("f", "")))
    out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
              % (g, m.get("sgb", ""), t.get("f", ""), t.get("g", ""),
                 m.get("n_mags", ""), m.get("rep_completeness", ""),
                 m.get("rep_contamination", ""),
                 r["msa_percent"], gpass, gcss, grrs, " | ".join(v["reasons"])))
out.close()

print()
print("  GUNC-failing among flagged: %d of %d" % (gunc_fail, len(flagged)))
print("  A genome flagged by BOTH GTDB-Tk and GUNC is a stronger exclusion")
print("  candidate than one flagged by either alone.")
print("  wrote %s" % OUT)

# ------------------------------------------------- are flagged genomes worse?
print()
print("=" * 74)
print("DO FLAGGED GENOMES DIFFER FROM THE REST?")
print("=" * 74)
wild = [r for r in man_rows if r["has_wild"].strip().lower() == "yes"]


def summ(vals, lab):
    vals = sorted(v for v in vals if v is not None)
    if not vals:
        print("    %s: none" % lab)
        return
    print("    %-22s n=%3d  min %6.2f  median %6.2f  max %6.2f"
          % (lab, len(vals), vals[0], vals[len(vals) // 2], vals[-1]))


for field, lab in (("rep_completeness", "completeness"),
                   ("rep_contamination", "contamination")):
    summ([tofloat(r[field]) for r in wild if r["representative"].strip() in flagged],
         "flagged " + lab)
    summ([tofloat(r[field]) for r in wild if r["representative"].strip() not in flagged],
         "unflagged " + lab)

print()
print("  flagged by family:")
for f, n in Counter(parse_tax(v["row"]["classification"]).get("f", "")
                    for v in flagged.values()).most_common():
    tot = sum(1 for r in new if parse_tax(r["classification"]).get("f", "") == f)
    print("     %-24s %3d of %4d (%.1f%%)" % (f, n, tot, pct(n, tot)))
print("  A family flagged well above the %.1f%% catalog rate is worth a look."
      % pct(len(flagged), len(new)))

# SENTINEL_END
