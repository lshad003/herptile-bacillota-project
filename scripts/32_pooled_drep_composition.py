# Dereplication clusters are resolved into one unit per arm per cluster.
# Source: /bigdata/stajichlab/lshad003/ruminococcaceae-agent/scripts/pooled_drep_composition.py
# Output: results/pooled_drep_cluster_composition.tsv
# POOLED_DREP_COMPOSITION_V2_20260804
import os
from collections import Counter, defaultdict

ROOT = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent"
WORK = os.path.join(ROOT, "work/pooled_drep_rum")
CDB = os.path.join(WORK, "drep_out/data_tables/Cdb.csv")
WDB = os.path.join(WORK, "drep_out/data_tables/Wdb.csv")
MAP = os.path.join(WORK, "genome_arms.tsv")
OUT = os.path.join(ROOT, "results/pooled_drep_cluster_composition.tsv")

ARMS_ORDER = ["herptile", "ehi_amphibian", "ehi", "youngblut", "gtdb_ref"]
ARMS = None  # set from genome_arms.tsv after it is read

arm_of = {}
gid_of = {}
with open(MAP) as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        arm_of[f[0]] = f[1]
        gid_of[f[0]] = f[2]
print("mapped genomes: %d  %s" % (len(arm_of), dict(Counter(arm_of.values()))))
seen = set(arm_of.values())
ARMS = [a for a in ARMS_ORDER if a in seen] + sorted(seen - set(ARMS_ORDER))
print("arms present: %s" % ARMS)

if not os.path.exists(CDB):
    raise SystemExit("no Cdb.csv at %s" % CDB)

clusters = defaultdict(list)
with open(CDB) as fh:
    hdr = fh.readline().rstrip("\n").replace('"', "").split(",")
    gi = hdr.index("genome")
    si = hdr.index("secondary_cluster")
    for line in fh:
        f = line.rstrip("\n").replace('"', "").split(",")
        if len(f) <= max(gi, si):
            continue
        clusters[f[si]].append(f[gi])

n_gen = sum(len(v) for v in clusters.values())
print()
print("=" * 74)
print("CLUSTERING")
print("=" * 74)
print("  genomes clustered : %d of %d staged" % (n_gen, len(arm_of)))
print("  dropped by -comp 70 -con 10 -l 50000: %d" % (len(arm_of) - n_gen))
print("  secondary clusters (cross-catalog SGBs): %d" % len(clusters))

surv = Counter(arm_of.get(g, "?") for v in clusters.values() for g in v)
print()
print("  survivors by arm:")
for a in ARMS:
    print("    %-10s %5d" % (a, surv.get(a, 0)))

# ---------------------------------------------- what is in each cluster
print()
print("=" * 74)
print("CLUSTER COMPOSITION: THE NOVELTY CHECK")
print("=" * 74)
print("  A cluster holding a herptile genome AND any reference or endotherm")
print("  genome means they are the same species (95%% ANI). R2 claims 100%%")
print("  species novelty for wild Ruminococcaceae, so that count must be 0.")

pat = Counter()
mixed = []
for cid, gs in clusters.items():
    a = set(arm_of.get(g, "?") for g in gs)
    key = "+".join(sorted(a))
    pat[key] += 1
    if "herptile" in a and len(a) > 1:
        mixed.append((cid, gs, a))

print()
print("  cluster arm-composition patterns:")
for k, v in pat.most_common():
    print("    %5d  %s" % (v, k))

print()
print("  clusters containing herptile AND another arm: %d" % len(mixed))
if mixed:
    print("  THIS CONTRADICTS THE 100%% NOVELTY CLAIM. Detail:")
    for cid, gs, a in mixed[:15]:
        print("    cluster %-10s arms=%s" % (cid, sorted(a)))
        for g in gs:
            print("        %-12s %s" % (arm_of.get(g, "?"), gid_of.get(g, g)))
else:
    print("  NONE. No wild herptile Ruminococcaceae genome shares a 95%% ANI")
    print("  species cluster with any endotherm MAG or GTDB reference.")
    print("  This is an independent confirmation of R2 at the genome level.")

# ------------------------------------------- herptile-only cluster count
herp_only = sum(1 for cid, gs in clusters.items()
                if set(arm_of.get(g, "?") for g in gs) == {"herptile"})
print()
print("  clusters composed only of herptile genomes: %d" % herp_only)

# ------------------------------------------- how redundant was ehi
print()
print("=" * 74)
print("HOW MUCH REDUNDANCY WAS IN THE UN-DEREPLICATED EHI ARM")
print("=" * 74)
for a in ARMS:
    genomes = sum(1 for v in clusters.values() for g in v if arm_of.get(g) == a)
    cl = len(set(cid for cid, v in clusters.items()
                 for g in v if arm_of.get(g) == a))
    if genomes:
        print("  %-10s %5d genomes in %5d clusters, %.2f genomes per cluster"
              % (a, genomes, cl, genomes / float(cl)))
print("  A ratio near 1 means the arm was already dereplicated. Higher means")
print("  that arm was contributing redundant near-identical genomes.")

with open(OUT, "w") as f:
    f.write("cluster\tn_genomes\tarms\therptile\tehi\tyoungblut\tgtdb_ref\n")
    for cid, gs in sorted(clusters.items()):
        c = Counter(arm_of.get(g, "?") for g in gs)
        f.write("%s\t%d\t%s\t%d\t%d\t%d\t%d\n"
                % (cid, len(gs), "+".join(sorted(set(c))),
                   c.get("herptile", 0), c.get("ehi", 0),
                   c.get("youngblut", 0), c.get("gtdb_ref", 0)))
print()
print("wrote %s" % OUT)

# SENTINEL_END
