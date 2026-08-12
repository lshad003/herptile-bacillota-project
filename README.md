# Herptile gut Bacillota_A: analysis record

Code and job scripts behind a comparative genomics study of Bacillota_A
metagenome-assembled genomes recovered from reptile and amphibian gut
metagenomes.

Paths inside the scripts point at the working repository, `ruminococcaceae-agent`,
and are left unchanged, so that a file here can be compared directly against the
version that was run. Genomes and intermediate data are not committed; they
remain on the UCR HPCC cluster.

## How the steps map to the manuscript

Steps are numbered by the order the analyses were run, not by manuscript
section. Filenames are numbered sequentially across the whole repository.

Section numbering is provisional and will change if the scope of the paper
changes.

| Repo step | Manuscript section |
|---|---|
| Step 1 | 3.1 catalog construction, and 3.3 assembly contiguity |
| Step 2 | 3.2 novelty against GTDB r220 |
| Step 3 | 3.4 genus composition against comparison catalogs |
| Step 4 | folded into 3.2 as the taxonomy testability result |
| Step 5 | 3.5 gene content |
| Step 6 | 3.6 biosynthetic gene clusters |

Step 4 was originally a full section on read profiling against genome
recovery. Only the count of genera whose names correspond between GTDB and
NCBI is retained in the manuscript; the scripts for the rest of that
analysis remain here because they were run.

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
| `scripts/25_make_step2_tables.py` | Supplementary tables S4 and S5 |

Output: `figures/Figure2_novelty_expansion.pdf`, `tables/TableS4`, `tables/TableS5`

### Step 3. Genus composition against comparison catalogs

Ruminococcaceae genome sets from two endotherm gut catalogs and a second,
independent amphibian catalog are classified against GTDB r220 on the same
taxonomy, then resolved into one unit per catalog per dereplication cluster so
that the four genus pools are counted in the same currency. Pool differences
are partitioned into turnover and nestedness, tested for sensitivity to the
evidence threshold, and checked against accumulation curves, since a genus
missing from a small catalog may reflect sampling rather than absence.

| File | Purpose |
|---|---|
| `scripts/26_stage_ehi_nonherptile.py` | EHI mammal arm staged |
| `jobs/27_gtdbtk_ehi_r220.sh` | EHI classification, identify and align |
| `jobs/28_gtdbtk_ehi_classify_resume.sh` | EHI classification, classify step |
| `scripts/29_stage_ehi_amphibian.py` | EHI newt arm staged |
| `jobs/30_gtdbtk_ehi_amphibian.sh` | EHI newt classification |
| `jobs/31_gtdbtk_youngblut.sh` | Youngblut classification |
| `scripts/32_pooled_drep_composition.py` | Clusters resolved to units per arm |
| `scripts/33_turnover_nestedness.py` | Turnover, nestedness, thresholds, accumulation |
| `scripts/34_four_catalog_genus_table.py` | Genus recovery blocks and matched Jaccard |
| `scripts/35_amphibian_genus_replication.py` | Recovery across the two amphibian catalogs |
| `scripts/36_figure3_cross_catalog.py` | Figure 3 |
| `scripts/37_make_step3_table.py` | Supplementary table S6 |

Output: `figures/Figure3_cross_catalog_overlap.pdf`, `tables/TableS6`

### Step 4. Read profiling against genome recovery

Genus-level read classifications from the wild samples are compared against the
genomes recovered from the same material. Because GTDB and NCBI genus names do
not correspond one to one, the comparison is restricted to genera whose names
map reciprocally between the two taxonomies, and genera failing that criterion
are reported as untestable rather than as undetected.

| File | Purpose |
|---|---|
| `scripts/38_gtdb_ncbi_genus_map.py` | GTDB genera mapped onto NCBI names |
| `scripts/39_gtdb_ncbi_reciprocal.py` | Reciprocal filter at 70 percent |
| `scripts/40_bracken_wild_bacillota.py` | Read profiles summarised by phylum and genus |
| `scripts/41_bracken_reciprocal_join.py` | Reads joined to genome recovery |
| `scripts/42_figure_reads_vs_genomes.py` | Reads versus genomes figure |
| `scripts/43_make_step4_table.py` | Supplementary table S7 |

Output: `figures/Figure_reads_vs_genomes.pdf`, `tables/TableS7`

### Step 5. Gene content where host and lineage can be separated

Gene content can be attributed to host origin only in genera where amphibian
and reference genomes are phylogenetically interleaved, so a joint tree is
built first and genera are tested for interleaving before any functional
comparison. Within the two qualifying genera, protein clusters are merged into
orthologous groups and carbohydrate-active enzyme families are annotated on an
independent layer. Both are modelled against genome completeness, since the
arms differ in it, and both are checked against label permutations. Per-protein
annotation rate is measured separately, because it sets the direction in which
any remaining bias runs.

| File | Purpose |
|---|---|
| `scripts/44_extract_ref_msa.py` | Reference alignments extracted and masked |
| `scripts/45_build_rep_tree_msa.py` | Query alignments assembled |
| `scripts/46_build_figure_tree_msa.py` | Tree input combined with outgroup |
| `jobs/47_figure_tree.sh` | Joint tree inferred |
| `scripts/48_add_ref_genus_to_tree_meta.py` | Genus assignments joined onto reference tips |
| `scripts/49_check_genus_interleaving.py` | Interleaving tested by Fitch parsimony |
| `scripts/50_stage_focal_genus_proteomes.py` | Focal proteomes staged |
| `jobs/51_mmseqs_focal_pangenomes.sh` | Protein clustering |
| `jobs/52_eggnog_focal.sh` | Orthologous group assignment |
| `scripts/53_merge_clusters_by_og.py` | Clusters merged by orthologous group |
| `jobs/54_checkm2_focal.sh` | Completeness and contamination |
| `scripts/55_run_happi_focal.R` | Clustering parameter comparison |
| `scripts/56_run_happi_merged.R` | Prevalence modelled against completeness |
| `jobs/57_happi_merged.sh` | Real fit and two permutations |
| `scripts/58_genus_filter_merged.py` | Direction required to agree within genera |
| `scripts/59_figure_gene_content.py` | Gene content figure |
| `scripts/60_stage_cazy_targets.py` | Proteomes staged for CAZy annotation |
| `jobs/61_cazy_hmmsearch.sh` | dbCAN search |
| `scripts/62_parse_cazy_focal.py` | Hits filtered and reduced to a family matrix |
| `scripts/63_run_happi_cazy.R` | CAZy prevalence modelled against completeness |
| `jobs/64_happi_cazy.sh` | Real fit and two permutations |
| `scripts/65_cazy_untested_families.py` | Rare families scanned across unassigned clades |
| `scripts/66_extract_pfam_models.py` | Individual Pfam models extracted |
| `jobs/67_annot_154.sh` | Pfam and KofamScan annotation |
| `scripts/68_parse_annot_154.py` | Annotation parsed into presence matrices |
| `scripts/69_annot_rate_by_arm.py` | Per-protein annotation rate compared between arms |
| `scripts/84_figure_cazy_tree_heatmap.py` | CAZy family heatmap |
| `scripts/85_make_step5_tables.py` | Supplementary tables S8 and S9 |

Output: `figures/Figure_gene_content.pdf`, `figures/Figure_cazy_tree_heatmap.pdf`, `tables/TableS8`, `tables/TableS9`

### Step 6. Biosynthetic gene clusters

Biosynthetic gene clusters are detected in all five arms and counted separately
as complete and as contig-edge, since a fragmented assembly splits clusters and
loses their boundaries. Cluster counts are then examined against assembly
contiguity within each arm before any between-arm statement is made, and product
class composition is compared only for classes whose proportion does not track
contiguity within arms.

| File | Purpose |
|---|---|
| `scripts/70_stage_bgc_input.py` | Amphibian genomes staged |
| `jobs/71_antismash_amphibian.sh` | antiSMASH, amphibian arm |
| `scripts/72_parse_antismash_amphibian.py` | Amphibian output parsed |
| `scripts/73_stage_bgc_refs.py` | Reference genomes staged |
| `jobs/74_antismash_refs.sh` | antiSMASH, reference arm |
| `scripts/75_parse_antismash_refs.py` | Reference output parsed |
| `scripts/76_stage_bgc_endotherm.py` | Comparison genomes staged |
| `jobs/77_antismash_endotherm.sh` | antiSMASH, comparison arms |
| `scripts/78_parse_antismash_endotherm.py` | Comparison output parsed |
| `scripts/79_assembly_quality_arms.py` | Contiguity measured from staged fastas |
| `scripts/80_assembly_quality_all_arms.py` | Comparison arms added |
| `scripts/81_bgc_density_and_carriage.py` | Density per megabase and carriage |
| `scripts/82_bgc_class_composition.py` | Class composition with the fragmentation gate |
| `scripts/83_figure_bgc_contiguity.py` | Figure |
| `scripts/86_make_step6_table.py` | Supplementary table S10 |

Output: `figures/Figure_bgc_contiguity.pdf`, `tables/TableS10`


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
