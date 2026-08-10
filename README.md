# Herptile gut Bacillota_A: analysis record

Code and job scripts behind a comparative genomics study of Bacillota_A
metagenome-assembled genomes recovered from reptile and amphibian gut
metagenomes.

Paths inside the scripts point at the working repository, `ruminococcaceae-agent`,
and are left unchanged, so that a file here can be compared directly against the
version that was run. Genomes and intermediate data are not committed; they
remain on the UCR HPCC cluster.

## Steps

### Step 1. Catalog construction and quality assessment

Bacillota_A MAGs recovered from herptile faecal metagenomes are dereplicated
into species-level genome bins, and each representative is assessed for
chimerism and for the rRNA and tRNA features required by the MIMAG standard.
Feature recovery is measured in a set of GTDB reference genomes as well, using
one parsing pipeline for both, so that the two are comparable.

| File | Purpose |
|---|---|
| `jobs/01_dereplicate_mags.sh` | Dereplication at 95% ANI |
| `scripts/02_build_sgb_manifest.py` | SGB manifest from the dereplication tables |
| `jobs/03_rrna_trna_catalog.sh` | barrnap and tRNAscan-SE, catalog arm |
| `jobs/05_rrna_trna_references.sh` | barrnap and tRNAscan-SE, reference arm |
| `scripts/06_parse_rrna_trna.py` | Feature parsing, both arms |
| `jobs/07_gunc_screen_representatives.sh` | Chimerism screening |
| `scripts/08_stage_gunc_input.py` | Protein calls for the above |
| `scripts/09_audit_gunc.py` | Chimerism summary |
| `scripts/10_figure1_catalog.py` | Figure 1 |
| `scripts/11_make_step1_tables.py` | Supplementary tables S1 to S3 |

Output: `figures/Figure1_catalog.pdf`, `tables/TableS1` to `TableS3`

### Step 2. Novelty against GTDB r220

Representatives of the wild species-level genome bins are classified against
GTDB r220 to establish how many fall outside existing species clusters, and
per-genus counts are compared against the species clusters GTDB already holds
for the same genus. The species-level result is then checked at genome level by
dereplicating the Ruminococcaceae representatives together with four comparison
sets in a single run, so that any shared species cluster would be recovered
directly rather than inferred from taxonomy.

| File | Purpose |
|---|---|
| `scripts/12_stage_wild_sgb_batchfile.py` | GTDB-Tk batch file for the wild representatives |
| `jobs/13_gtdbtk_wild_sgb.sh` | Classification against GTDB r220 |
| `scripts/14_verify_wild_classify.py` | Assignments compared against the SGB manifest |
| `scripts/15_audit_wild_classify_flags.py` | GTDB-Tk quality flags against contamination and GUNC |
| `scripts/16_novelty_proportions.py` | Novelty proportions and per-genus expansion |
| `scripts/17_sampling_effort_all_genera.py` | Expansion against GTDB species-cluster counts |
| `scripts/18_genus_r220_to_r226_map.py` | Genus stability from r220 to r226 |
| `scripts/19_stage_rum_drep.py` | Dereplication input set, first arm |
| `scripts/20_stage_gtdb_refs_rum.py` | Dereplication input set, GTDB references |
| `scripts/21_stage_amph_into_drep.py` | Dereplication input set, amphibian arms |
| `jobs/22_pooled_drep_rum.sh` | Joint dereplication of five arms at 95% ANI |
| `scripts/23_pooled_drep_multiarm_audit.py` | Cluster composition by arm |
| `scripts/24_figure2_novelty.py` | Figure 2 and the genus-expansion permutation |

Output: `figures/Figure2_novelty_expansion.pdf`

## Software

GTDB-Tk 2.4.1 (GTDB r220), dRep 3.5.0, CheckM1, CheckM2 1.1.0, GUNC 1.0.6,
MMseqs2, Prodigal V2.6.3, FastTree 2.1.11, HMMER 3.4, Pfam 35.0 and 38.2,
KofamScan 1.3.0 (KEGG 97.0), dbCAN v13.0, eggNOG-mapper, antiSMASH 7.1.0,
barrnap 0.9, tRNAscan-SE 2.0.12, happi (R).

## Repository layout

    jobs/      SLURM submission scripts
    scripts/   analysis code
    figures/   final manuscript figures
    tables/    supplementary tables
    config/    input lists and batch files
    docs/      extended notes
