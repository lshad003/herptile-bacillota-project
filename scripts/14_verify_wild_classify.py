# GTDB-Tk assignments for the wild representatives are compared against the SGB manifest.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/verify_wild_classify.py
# Output: results/wild_classify_verification.tsv
# VERIFY_WILD_CLASSIFY_V1_20260803
import os
from collections import Counter

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
SUM = os.path.join(ROOT, "results/gtdbtk_wild_sgb_r220/gtdbtk.bac120.summary.tsv")
MANIFEST = os.path.join(ROOT, "data/sgb_manifest.tsv")
OUT = os.path.join(ROOT, "results/wild_classify_verification.tsv")


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


def norm(g):
    g = g.strip()
    return "" if g.upper() in ("UNASSIGNED", "", "NA", "N/A") else g


def pct(a, b):
    return (100.0 * a / b) if b else 0.0


new = read_tsv(SUM)
man = {r["representative"].strip(): r for r in read_tsv(MANIFEST)
       if r["has_wild"].strip().lower() == "yes"}

print("=" * 74)
print("COVERAGE")
print("=" * 74)
print("  summary rows: %d | wild SGBs in manifest: %d" % (len(new), len(man)))
ids = {r["user_genome"].strip() for r in new}
print("  matched: %d | in summary not manifest: %d | in manifest not summary: %d"
      % (len(ids & set(man)), len(ids - set(man)), len(set(man) - ids)))

# ------------------------------------------------------------- warnings
print()
print("=" * 74)
print("THE 28 WARNINGS")
print("=" * 74)
warn = [r for r in new if r["warnings"].strip() and r["warnings"].strip() != "N/A"]
print("  genomes with a warning: %d" % len(warn))
for w, n in Counter(r["warnings"].strip() for r in warn).most_common():
    print("    %4d  %s" % (n, w[:96]))
if warn:
    print("  families affected: %s"
          % dict(Counter(parse_tax(r["classification"]).get("f", "") for r in warn)))

# ------------------------------------------------- classification method
print()
print("=" * 74)
print("CLASSIFICATION METHOD")
print("=" * 74)
for m, n in Counter(r["classification_method"].strip() for r in new).most_common():
    print("  %4d (%5.1f%%)  %s" % (n, pct(n, len(new)), m[:80]))

msa = [tofloat(r["msa_percent"]) for r in new]
msa = sorted(v for v in msa if v is not None)
if msa:
    print("  msa_percent: min %.1f  median %.1f  max %.1f  below 50%%: %d"
          % (msa[0], msa[len(msa) // 2], msa[-1], sum(1 for v in msa if v < 50)))

# --------------------------------------------- agreement with the manifest
print()
print("=" * 74)
print("AGREEMENT WITH data/sgb_manifest.tsv")
print("=" * 74)
print("  The manifest genus came from a run whose summary is gone. This is the")
print("  first time those assignments can be checked against a live file.")
fam_d = []
gen_d = []
sp_d = []
for r in new:
    g = r["user_genome"].strip()
    if g not in man:
        continue
    t = parse_tax(r["classification"])
    m = man[g]
    if t.get("f", "") != m["family"].strip():
        fam_d.append((g, m["family"].strip(), t.get("f", "")))
    if norm(t.get("g", "")) != norm(m["genus"]):
        gen_d.append((g, norm(m["genus"]), norm(t.get("g", ""))))
    if t.get("s", "") != m["species"].strip():
        sp_d.append((g, m["species"].strip(), t.get("s", "")))
print("  family disagreements : %d" % len(fam_d))
for g, a, b in fam_d[:10]:
    print("      %-36s manifest %-22s new %s" % (g, a or "<blank>", b or "<blank>"))
print("  genus disagreements  : %d" % len(gen_d))
for g, a, b in gen_d[:15]:
    print("      %-36s manifest %-22s new %s" % (g, a or "<blank>", b or "<blank>"))
print("  species disagreements: %d" % len(sp_d))
for g, a, b in sp_d[:10]:
    print("      %-36s manifest %-22s new %s" % (g, a or "<blank>", b or "<blank>"))

# ------------------------------------------------ novelty, ANI-verified
print()
print("=" * 74)
print("NOVELTY, NOW VERIFIABLE AGAINST ANI AND RADIUS")
print("=" * 74)
out = open(OUT, "w")
out.write("genome\tfamily\tmanifest_genus\tnew_genus\tspecies\tani\taf\tradius\tverdict\n")
blank = 0
no_ani = 0
below = 0
contra = 0
for r in new:
    g = r["user_genome"].strip()
    t = parse_tax(r["classification"])
    sp = t.get("s", "")
    ani = tofloat(r["closest_genome_ani"])
    af = tofloat(r["closest_genome_af"])
    rad = tofloat(r["closest_genome_reference_radius"])
    if sp:
        continue
    blank += 1
    if ani is None:
        no_ani += 1
        v = "no_ani_consistent"
    elif rad is not None and ani >= rad and (af is not None and af >= 0.5):
        contra += 1
        v = "CONTRADICTION"
    else:
        below += 1
        v = "below_radius_consistent"
    out.write("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\n"
              % (g, t.get("f", ""), norm(man.get(g, {}).get("genus", "")),
                 norm(t.get("g", "")), r["closest_genome_ani"],
                 r["closest_genome_af"],
                 r["closest_genome_reference_radius"], v))
out.close()
print("  SGBs with no species assignment: %d of %d (%.1f%%)"
      % (blank, len(new), pct(blank, len(new))))
print("    no ANI reported, consistent  : %d" % no_ani)
print("    ANI below radius, consistent : %d" % below)
print("    CONTRADICTIONS               : %d" % contra)
print("  wrote %s" % OUT)

# ------------------------------------------------ R2 numbers recomputed
print()
print("=" * 74)
print("R2 NUMBERS, RECOMPUTED FROM THE NEW SUMMARY")
print("=" * 74)
rum = [r for r in new if parse_tax(r["classification"]).get("f", "") == "Ruminococcaceae"]
for lab, sub, e_n, e_sp, e_gen in (
        ("all wild SGBs", new, 718, 715, 215),
        ("wild Ruminococcaceae", rum, 220, 220, 25)):
    nsp = sum(1 for r in sub if not parse_tax(r["classification"]).get("s", ""))
    ngen = sum(1 for r in sub if not parse_tax(r["classification"]).get("g", ""))
    nfam = sum(1 for r in sub if not parse_tax(r["classification"]).get("f", ""))
    print("  %s: n=%d (expect %d)" % (lab, len(sub), e_n))
    print("     no species: %d (expect %d)  %s"
          % (nsp, e_sp, "OK" if nsp == e_sp else "CHECK"))
    print("     no genus  : %d (expect %d)  %s"
          % (ngen, e_gen, "OK" if ngen == e_gen else "CHECK"))
    print("     no family : %d" % nfam)

print()
print("  genus counts, new summary, wild Ruminococcaceae:")
for g, n in Counter(norm(parse_tax(r["classification"]).get("g", ""))
                    for r in rum).most_common(12):
    print("     %-28s %4d" % (g if g else "<no genus>", n))

# SENTINEL_END
