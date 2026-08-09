# Herptile gut Bacillota_A: analysis record

Code and job scripts behind a comparative genomics study of Bacillota_A
metagenome-assembled genomes recovered from reptile and amphibian gut
metagenomes.

Scripts are copied from the working repository, `ruminococcaceae-agent`, and
only those that produced reported results are included here. Each file records
its original path, what it reads, and what it writes. Genomes and intermediate
data are not committed; they remain on the UCR HPCC cluster.

Paths inside the scripts point at the working repository and are left
unchanged, so that a file here can be compared directly against the version
that was run.

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
