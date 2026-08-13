# Herptile gut Bacillota_A: analysis record

This repository documents the analyses used to characterize Bacillota_A
genomes recovered from reptile and amphibian gut metagenomes, determine how
much of that diversity is represented in existing genome databases, compare
Ruminococcaceae composition across independent genome catalogs, and examine
gene-content differences both where host origin can be separated from
bacterial lineage and where the two remain confounded, alongside
biosynthetic feature comparisons.

Each step below summarizes the analysis, its main result, the scripts used,
and the resulting figures or supplementary tables.

## Steps

### Step 1. Catalog construction and quality assessment

Bacillota_A MAGs recovered from herptile faecal metagenomes are dereplicated
into species-level genome bins, and each representative is assessed for
chimerism and for the rRNA and tRNA features required by the MIMAG standard.
Feature recovery is measured in a set of GTDB reference genomes as well, using
one parsing pipeline for both, so that the two are comparable.

**Main result.** The catalog contains 2,229 Bacillota_A MAGs, dereplicated at 95% ANI to 1,171 SGB representatives, including 718 with at least one wild-derived MAG. Genome-quality screens supported the integrity of the catalog. Comparison with GTDB Ruminococcaceae references showed that short genomic features are recovered at comparable rates in both sets, while recovery of longer features tracks assembly contiguity, which sets the expectation for the long-feature analyses later in the study.

| File | Purpose |
|---|---|
| `jobs/01_dereplicate_mags.sh` | Dereplication at 95% ANI |
| `scripts/02_build_sgb_manifest.py` | SGB manifest from the dereplication tables |
| `jobs/03_rrna_trna_catalog.sh` | barrnap and tRNAscan-SE, catalog set |
| `jobs/05_rrna_trna_references.sh` | barrnap and tRNAscan-SE, reference set |
| `scripts/06_parse_rrna_trna.py` | Feature parsing, both sets |
| `jobs/07_gunc_screen_representatives.sh` | Chimerism screening |
| `scripts/08_stage_gunc_input.py` | Protein calls for the above |
| `scripts/09_audit_gunc.py` | Chimerism summary |
| `scripts/10_figure1_catalog.py` | Figure 1 |
| `scripts/11_make_step1_tables.py` | Supplementary tables S1 to S3 |
| `scripts/87_figure_feature_length_ratio.py` | Feature recovery against feature length |

Output: `figures/Figure1_catalog.pdf`, `figures/Figure_feature_length_ratio.pdf`, `tables/TableS1` to `TableS3`

### Step 2. Novelty against GTDB r220

Representatives of the wild species-level genome bins are classified against
GTDB r220 to establish how many fall outside existing species clusters, and
per-genus counts are compared against the species clusters GTDB already holds
for the same genus. The species-level result is then checked at genome level by
dereplicating the Ruminococcaceae representatives together with four comparison
sets in a single run, so that any shared species cluster would be recovered
directly rather than inferred from taxonomy.

**Main result.** Of the 718 wild SGB representatives, 715 (99.6%) had no match to an existing GTDB r220 species cluster. Joint dereplication at 95% ANI with four external genome sets produced no cluster containing both a wild-catalog genome and a genome from any other set, while 57 clusters contained genomes from more than one external set, so cross-set clustering did occur in the same run. Within Ruminococcaceae, expansion was concentrated in particular genera: eight held at least twice as many wild SGBs as GTDB r220 species clusters, whereas no more than five genera reached that level in any of 9,999 permutations that redistributed the same SGBs in proportion to existing GTDB representation (p = 0.0001). The twofold threshold was chosen after inspecting the data, so this is exploratory. Novelty also extended above species level: 215 of 718 wild representatives (29.9%), including 25 of the 220 wild Ruminococcaceae, had no GTDB genus assignment, and relative evolutionary divergence supported these as deeper phylogenetic placements rather than classification failures. No names are proposed, since the genomes placing deepest are also the least contiguous.

| File | Purpose |
|---|---|
| `scripts/12_stage_wild_sgb_batchfile.py` | GTDB-Tk batch file for the wild representatives |
| `jobs/13_gtdbtk_wild_sgb.sh` | Classification against GTDB r220 |
| `scripts/14_verify_wild_classify.py` | Assignments compared against the SGB manifest |
| `scripts/15_audit_wild_classify_flags.py` | GTDB-Tk quality flags against contamination and GUNC |
| `scripts/16_novelty_proportions.py` | Novelty proportions and per-genus expansion |
| `scripts/17_sampling_effort_all_genera.py` | Expansion against GTDB species-cluster counts |
| `scripts/18_genus_r220_to_r226_map.py` | Genus stability from r220 to r226 |
| `scripts/19_stage_rum_drep.py` | Dereplication input, first set |
| `scripts/20_stage_gtdb_refs_rum.py` | Dereplication input set, GTDB references |
| `scripts/21_stage_amph_into_drep.py` | Dereplication input, amphibian sets |
| `jobs/22_pooled_drep_rum.sh` | Joint dereplication of five sets at 95% ANI |
| `scripts/23_pooled_drep_multiarm_audit.py` | Cluster composition by set |
| `scripts/24_figure2_novelty.py` | Figure 2 and the genus-expansion permutation |
| `scripts/25_make_step2_tables.py` | Supplementary tables S4 and S5 |
| `scripts/88_unassigned_clade_coherence.py` | Genus-unassigned tips grouped into clades |
| `scripts/89_red_genus_check.py` | Relative evolutionary divergence of unassigned genomes |
| `scripts/90_figure_family_tree.py` | Family tree with divergence and assignment tracks |

Output: `figures/Figure2_novelty_expansion.pdf`, `figures/Figure_family_tree_red.pdf`, `tables/TableS4`, `tables/TableS5`

**
### Step 3. Genus composition against comparison catalogs

Ruminococcaceae genome sets from the EHI and Youngblut catalogs and a second
amphibian catalog are classified against GTDB r220 on the same
taxonomy, then resolved into one unit per catalog per dereplication cluster so
that the four genus pools are counted in the same currency. Pool differences
are partitioned into turnover and nestedness, tested for sensitivity to the
evidence threshold, and checked against accumulation curves, since a genus
missing from a small catalog may reflect sampling rather than absence.

**Main result.** Ruminococcaceae genus composition in the wild catalog differs from the EHI comparison catalog by replacement rather than nested loss. Both sets hold 23 genera and share four, giving Sørensen dissimilarity 0.826 with no nestedness component. The EHI newt catalog, which samples different animals but shares the EHI pipeline, is closer to the wild catalog than to the EHI set, with all 14 EHI newt genera represented among the 23 wild genera. No set has saturated, so absences in the smaller catalogs should be interpreted cautiously.

| File | Purpose |
|---|---|
| `scripts/26_stage_ehi_nonherptile.py` | EHI mammal set staged |
| `jobs/27_gtdbtk_ehi_r220.sh` | EHI classification, identify and align |
| `jobs/28_gtdbtk_ehi_classify_resume.sh` | EHI classification, classify step |
| `scripts/29_stage_ehi_amphibian.py` | EHI newt set staged |
| `jobs/30_gtdbtk_ehi_amphibian.sh` | EHI newt classification |
| `jobs/31_gtdbtk_youngblut.sh` | Youngblut classification |
| `scripts/32_pooled_drep_composition.py` | Clusters resolved to units per set |
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

**Main result.** Comparison between genome recovery and NCBI-based read profiling is limited by taxonomy rather than by the number of genomes recovered. Of 123 genera in the wild catalog, only 27 have reciprocal GTDB to NCBI genus mappings permitting an unambiguous comparison; 39 of the remainder carry placeholder names with no NCBI equivalent on any reference genome. Among the 27, 14 were classified in reads and 13 were recovered as genomes without any genus-level read assignment. This step measures an interoperability limit between two taxonomies rather than providing a complete abundance comparison.

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
sets differ in it, and both are checked against label permutations. Per-protein
annotation rate is measured separately, because it sets the direction in which
any remaining bias runs.

**Main result.** Gene content comparisons were restricted to the two genera in which amphibian and reference genomes are phylogenetically interleaved, so that host origin is not confounded with lineage. Of 3,345 orthologous groups constructed, 3,278 were fitted and 618 differed in detectable prevalence after correction for genome completeness. Two label permutations returned zero of 3,278. Carbohydrate-active enzyme families were tested independently on the same genomes, with 22 of 196 testable families differing and both permutations returning zero.

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
| `scripts/69_annot_rate_by_arm.py` | Per-protein annotation rate compared between sets |
| `scripts/84_figure_cazy_tree_heatmap.py` | CAZy family heatmap |
| `scripts/85_make_step5_tables.py` | Supplementary tables S8 and S9 |

Output: `figures/Figure_gene_content.pdf`, `figures/Figure_cazy_tree_heatmap.pdf`, `tables/TableS8`, `tables/TableS9`

### Step 6. Gene content where host and lineage cannot be separated
In Angelakisella, all 31 amphibian genomes (21 wild-catalog genomes and 10
EHI newt genomes) form a single clade within the GTDB-defined genus,
separate from 26 non-amphibian genomes. Host origin and bacterial lineage
therefore cannot be separated in this genus. Rather than testing for a
host-origin association, this analysis asks whether the replicated
amphibian-associated lineage carries a distinct functional repertoire,
described throughout as clade-associated. Because sequence-identity
clustering can split divergent orthologs, apparent clade-specific MMseqs
clusters were first audited by homology search against the reference
proteomes, using clusters shared between groups as a positive control. Gene
content was then represented by eggNOG orthologous groups and tested with
the same completeness-aware framework and label-permutation controls used
for the interleaved genera.

Main result. Of 67 MMseqs protein clusters that were near-fixed in the
amphibian clade and near-absent from the other Angelakisella genomes, 66
had detectable reference homologs; their median best-hit identity was 56%,
closely matching the measured mean amphibian-to-reference AAI of 55.8%
(range 52.1 to 59.1% across 100 genome pairs). Across these genomes, MMseqs
clustering produced 3.1-fold more units than ortholog-level representation,
showing that much of the apparent clade-specific protein content reflected
sequence divergence rather than gene gain or loss. After orthology
correction and completeness-aware testing, 866 of 3,023 fitted orthologous
groups differed between the amphibian clade and other Angelakisella genomes
(504 amphibian-enriched and 362 reference-enriched), while both
permuted-label controls produced zero significant groups. Eight orthologous
groups were near-fixed in the amphibian clade and near-absent elsewhere,
including an adjacent kdpAB potassium-transport locus and an FkbH-like
family. Reference-enriched groups also showed coordinated differences in
biosynthetic functions, including chorismate- and pyrimidine-synthesis
genes.

File	Purpose
`scripts/91_angelakisella_matrices.py`	Neighborhood genome set and presence matrices in two currencies
`scripts/92_stage_diamond_queries.py`	Cluster member proteins staged for homology search
`jobs/93_angelakisella_diamond.sh`	Homology search against reference proteomes
`scripts/94_classify_diamond_hits.py`	Clusters classified with shared-cluster positive control
`scripts/95_nohit_family_members.py`	Membership of the family without reference homologs
`jobs/96_nohit_family_tblastn.sh`	Absence closed by translated search against reference assemblies
`scripts/97_nohit_family_hmmscan.sh`	Domain identification of the family
`scripts/98_stage_emapper_gap.py`	Unannotated genomes staged
`jobs/99_angelakisella_emapper.sh`	eggNOG annotation, first set
`scripts/100_stage_emapper_gap2.py`	Additional genomes staged
`jobs/101_angelakisella_emapper2.sh`	eggNOG annotation, second set
`scripts/102_happi_metadata.py`	Completeness joined onto test genomes
`scripts/103_run_happi_angelakisella.R`	Prevalence modelled against completeness
`jobs/104_angelakisella_happi.sh`	Real fit and two permutations
`scripts/105_block_privacy.py`	Private units per clade in both currencies
`scripts/106_kdp_verify.py`	Kdp locus verified from annotations and adjacency
`scripts/107_refenriched_check.py`	Reference-enriched groups annotated from member proteins
`jobs/108_angelakisella_aai.sh`	Amphibian-to-reference AAI over 100 genome pairs
`scripts/109_figure_angelakisella.py`	Neighborhood tree and two-currency heatmap
Output: `figures/Figure_angelakisella_neighborhood_heatmap.pdf`


### Step 7. Biosynthetic gene clusters

Biosynthetic gene clusters are detected in all five sets and counted separately
as complete and as contig-edge, since a fragmented assembly splits clusters and
loses their boundaries. Cluster counts are then examined against assembly
contiguity within each set before any between-set statement is made, and product
class composition is compared only for classes whose proportion does not track
contiguity within sets.

**Main result.** Total antiSMASH region recovery is similar across the five Ruminococcaceae genome sets, whereas recovery of complete clusters depends on assembly contiguity. Complete-region counts track contig N50 within every set (Spearman rho 0.54 to 0.73) while total-region counts do not (rho 0.01 to 0.25), and at matched N50 of 40 to 80 kb the wild amphibian, EHI newt and GTDB sets give 0.67, 0.65 and 0.68 complete regions per genome. After the fragmentation gate, RRE-containing clusters are the only retained class showing the same direction in both amphibian sets. These comparisons are descriptive and do not establish biosynthetic novelty.

| File | Purpose |
|---|---|
| `scripts/70_stage_bgc_input.py` | Amphibian genomes staged |
| `jobs/71_antismash_amphibian.sh` | antiSMASH, amphibian set |
| `scripts/72_parse_antismash_amphibian.py` | Amphibian output parsed |
| `scripts/73_stage_bgc_refs.py` | Reference genomes staged |
| `jobs/74_antismash_refs.sh` | antiSMASH, reference set |
| `scripts/75_parse_antismash_refs.py` | Reference output parsed |
| `scripts/76_stage_bgc_endotherm.py` | Comparison genomes staged |
| `jobs/77_antismash_endotherm.sh` | antiSMASH, comparison sets |
| `scripts/78_parse_antismash_endotherm.py` | Comparison output parsed |
| `scripts/79_assembly_quality_arms.py` | Contiguity measured from staged fastas |
| `scripts/80_assembly_quality_all_arms.py` | Comparison sets added |
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
