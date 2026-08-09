#!/usr/bin/env python3
# Figure 1: sampling, genome quality, dereplication by family, chimerism.
#
# Source: ruminococcaceae-agent/scripts/fig1_catalog.py (V2)
# Reads:  data/herptile_bacillota_A_HQ_manifest_with_source.tsv
#         data/sgb_manifest.tsv
#         results/gunc_audit_by_sgb.tsv
# Writes: results/figures/Figure1_catalog.pdf and .png
#
# Panel a counts MAGs and is broken down by collection source. Panels c and d
# count SGBs. The two units differ and are labelled accordingly, since MAG
# counts and SGB counts diverge most where a single organism was recovered
# repeatedly from the same animals.
